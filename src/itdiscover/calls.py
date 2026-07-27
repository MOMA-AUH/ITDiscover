"""Build reportable ITD calls from aligned reads."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field, replace

from .alleles import (
    CanonicalInsertionAllele,
    canonicalize_insertion,
    canonicalize_insertion_allele,
)
from .coverage import (
    observed_mutant_fragment_fraction,
    spans_insertion_site,
)
from .insertions import (
    Alignment,
    Insertion,
    InsertionEvidenceFilter,
    extract_insertions,
)
from .itds import ITD, classify_exact_itd, classify_fuzzy_itd


@dataclass(frozen=True)
class ITDCall:
    """An ITD with mutually exclusive candidate-specific fragment evidence."""

    itd: ITD
    canonical_allele: CanonicalInsertionAllele
    mutant_fragment_count: int
    wild_type_fragment_count: int
    status: str = "PASS"
    filter_reasons: tuple[str, ...] = ()
    r1_mutant_count: int = field(default=0, compare=False)
    r2_mutant_count: int = field(default=0, compare=False)
    r1_opportunity_count: int = field(default=0, compare=False)
    r2_opportunity_count: int = field(default=0, compare=False)
    concordant_fragment_count: int = field(default=0, compare=False)
    single_mate_fragment_count: int = field(default=0, compare=False)
    conflicting_fragment_count: int = field(default=0, compare=False)
    unresolved_fragment_count: int = field(default=0, compare=False)
    not_informative_fragment_count: int = field(default=0, compare=False)
    consolidated_members: tuple["ConsolidatedAlleleMember", ...] = field(
        default=(),
        compare=False,
    )

    def __post_init__(self) -> None:
        counts = (
            self.mutant_fragment_count,
            self.wild_type_fragment_count,
            self.conflicting_fragment_count,
            self.unresolved_fragment_count,
            self.not_informative_fragment_count,
            self.r1_mutant_count,
            self.r2_mutant_count,
            self.r1_opportunity_count,
            self.r2_opportunity_count,
        )
        if any(count < 0 for count in counts):
            raise ValueError("fragment evidence counts must not be negative")
        if self.r1_mutant_count > self.r1_opportunity_count:
            raise ValueError(
                "r1_mutant_count must not exceed r1_opportunity_count"
            )
        if self.r2_mutant_count > self.r2_opportunity_count:
            raise ValueError(
                "r2_mutant_count must not exceed r2_opportunity_count"
            )

    @property
    def informative_fragment_count(self) -> int:
        """Return mutant plus high-quality wild-type fragments."""
        return self.mutant_fragment_count + self.wild_type_fragment_count

    @property
    def observed_mutant_fragment_fraction(self) -> float:
        """Return mutant fragments divided by informative fragments."""
        return observed_mutant_fragment_fraction(
            self.mutant_fragment_count,
            self.informative_fragment_count,
        )

    @property
    def r1_mutant_fraction(self) -> float | None:
        """Return the mutant fraction among R1 junction opportunities."""
        if not self.r1_opportunity_count:
            return None
        return self.r1_mutant_count / self.r1_opportunity_count

    @property
    def r2_mutant_fraction(self) -> float | None:
        """Return the mutant fraction among R2 junction opportunities."""
        if not self.r2_opportunity_count:
            return None
        return self.r2_mutant_count / self.r2_opportunity_count

    @property
    def passes_filters(self) -> bool:
        """Return whether the call passes the configured thresholds."""
        return self.status == "PASS"

    @property
    def consolidated_minor_fragment_count(self) -> int:
        """Return raw passing-fragment support absorbed from minor alleles."""
        return sum(member.fragment_count for member in self.consolidated_members)


@dataclass(frozen=True)
class ITDFilter:
    """Thresholds used to label ITD calls."""

    min_mutant_fragment_count: int = 3
    min_informative_fragment_count: int = 10
    min_observed_mutant_fragment_fraction: float = 0.01
    max_directional_mutant_fraction_share: float = 0.90
    min_directional_opportunities: int = 5

    def __post_init__(self) -> None:
        if self.min_mutant_fragment_count < 1:
            raise ValueError("min_mutant_fragment_count must be at least 1")
        if self.min_informative_fragment_count < 0:
            raise ValueError("min_informative_fragment_count must not be negative")
        if not 0 <= self.min_observed_mutant_fragment_fraction <= 1:
            raise ValueError(
                "min_observed_mutant_fragment_fraction must be between 0 and 1"
            )
        if not 0.5 <= self.max_directional_mutant_fraction_share <= 1:
            raise ValueError(
                "max_directional_mutant_fraction_share must be between 0.5 and 1"
            )
        if self.min_directional_opportunities < 1:
            raise ValueError("min_directional_opportunities must be at least 1")


@dataclass(frozen=True)
class ITDConsolidationSettings:
    """Safeguards for merging weak observations into a dominant ITD allele."""

    enabled: bool = False
    max_allele_mismatch_rate: float = 0.125
    max_breakpoint_shift_rate: float = 1.0
    max_minor_to_anchor_support_ratio: float = 0.05
    min_anchor_fragment_count: int = 3

    def __post_init__(self) -> None:
        if not 0 <= self.max_allele_mismatch_rate <= 1:
            raise ValueError(
                "max_allele_mismatch_rate must be between 0 and 1"
            )
        if not 0 <= self.max_breakpoint_shift_rate <= 1:
            raise ValueError(
                "max_breakpoint_shift_rate must be between 0 and 1"
            )
        if not 0 <= self.max_minor_to_anchor_support_ratio <= 1:
            raise ValueError(
                "max_minor_to_anchor_support_ratio must be between 0 and 1"
            )
        if self.min_anchor_fragment_count < 1:
            raise ValueError("min_anchor_fragment_count must be at least 1")


@dataclass(frozen=True)
class ConsolidatedAlleleMember:
    """One minor observed allele absorbed into a dominant call."""

    allele: CanonicalInsertionAllele
    fragment_count: int
    allele_mismatches: int
    allele_mismatch_rate: float
    breakpoint_shift: int
    breakpoint_shift_rate: float
    reason: str


@dataclass(frozen=True)
class UniqueSupportRepresentative:
    """One representative alignment for a unique local ITD support pattern."""

    itd: ITD
    signature: str
    alignment: Alignment
    canonical_allele: CanonicalInsertionAllele
    support_count: int
    exact_support_count: int
    fuzzy_only_support_count: int = 0
    fuzzy_example_sequence: str | None = None
    mismatches: int = 0
    insert_sequence_supports: tuple["InsertSequenceSupport", ...] = ()


@dataclass(frozen=True)
class InsertSequenceSupport:
    """Fragment support for one observed inserted sequence."""

    sequence: str
    support_count: int
    mismatches: int


@dataclass(frozen=True)
class FragmentConsensusSupport:
    """Mutually exclusive fragment evidence states for one candidate."""

    mutant_fragment_ids: frozenset[str] = frozenset()
    wild_type_fragment_ids: frozenset[str] = frozenset()
    concordant_fragment_ids: frozenset[str] = frozenset()
    single_mate_fragment_ids: frozenset[str] = frozenset()
    conflicting_fragment_ids: frozenset[str] = frozenset()
    unresolved_fragment_ids: frozenset[str] = frozenset()
    not_informative_fragment_ids: frozenset[str] = frozenset()
    r1_mutant_fragment_ids: frozenset[str] = frozenset()
    r2_mutant_fragment_ids: frozenset[str] = frozenset()
    r1_opportunity_fragment_ids: frozenset[str] = frozenset()
    r2_opportunity_fragment_ids: frozenset[str] = frozenset()


@dataclass(frozen=True)
class _CandidateObservation:
    itd: ITD
    alignment: Alignment
    observed_allele: CanonicalInsertionAllele
    canonical_allele: CanonicalInsertionAllele
    passes_evidence: bool
    consolidated_from: CanonicalInsertionAllele | None = None


# Exact observed ALT identity is independent of an aligner's arbitrary gap
# placement and of the copied-segment/spacer annotation inferred from that gap.
ITDCallKey = CanonicalInsertionAllele
SupportRepresentativeMap = dict[ITDCallKey, dict[str, UniqueSupportRepresentative]]


def call_exact_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
    min_copied_segment_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
    consolidation: ITDConsolidationSettings = ITDConsolidationSettings(),
) -> list[ITDCall]:
    """Call exact-match ITDs and attach fragment counts and observed fraction."""
    calls, _ = call_exact_itds_with_representatives(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        min_copied_segment_length=min_copied_segment_length,
        require_in_frame=require_in_frame,
        filters=filters,
        evidence_filter=evidence_filter,
        consolidation=consolidation,
    )
    return calls


def call_fuzzy_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_copy_mismatch_rate: float,
    min_insert_length: int = 6,
    min_copied_segment_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
    consolidation: ITDConsolidationSettings = ITDConsolidationSettings(),
) -> list[ITDCall]:
    """Call fuzzy-match ITDs and attach fragment counts and observed fraction."""
    if not 0 <= max_copy_mismatch_rate <= 1:
        raise ValueError("max_copy_mismatch_rate must be between 0 and 1")
    calls, _ = call_fuzzy_itds_with_representatives(
        alignments,
        reference,
        max_copy_mismatch_rate=max_copy_mismatch_rate,
        min_insert_length=min_insert_length,
        min_copied_segment_length=min_copied_segment_length,
        require_in_frame=require_in_frame,
        filters=filters,
        evidence_filter=evidence_filter,
        consolidation=consolidation,
    )
    return calls


def call_fuzzy_itds_with_representatives(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_copy_mismatch_rate: float,
    min_insert_length: int = 6,
    min_copied_segment_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
    consolidation: ITDConsolidationSettings = ITDConsolidationSettings(),
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call fuzzy-match ITDs and retain one alignment per unique support pattern."""
    if not 0 <= max_copy_mismatch_rate <= 1:
        raise ValueError("max_copy_mismatch_rate must be between 0 and 1")
    min_copied_segment_length = _resolved_min_copied_segment_length(
        min_insert_length,
        min_copied_segment_length,
    )
    alignments = list(alignments)
    (
        grouped_itds,
        representative_map,
        consensus_by_key,
        consolidated_members_by_key,
    ) = _collect_fuzzy_itd_support(
        alignments,
        reference,
        max_copy_mismatch_rate=max_copy_mismatch_rate,
        min_insert_length=min_insert_length,
        min_copied_segment_length=min_copied_segment_length,
        require_in_frame=require_in_frame,
        evidence_filter=evidence_filter,
        consolidation=consolidation,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for key, itds in grouped_itds.items():
        representative = _representative_itd(itds, key)
        consensus = consensus_by_key[key]
        mutant_fragment_count = len(consensus.mutant_fragment_ids)
        wild_type_fragment_count = len(consensus.wild_type_fragment_ids)
        informative_fragment_count = (
            mutant_fragment_count + wild_type_fragment_count
        )
        r1_mutant_count = len(consensus.r1_mutant_fragment_ids)
        r2_mutant_count = len(consensus.r2_mutant_fragment_ids)
        r1_opportunity_count = len(
            consensus.r1_opportunity_fragment_ids
        )
        r2_opportunity_count = len(
            consensus.r2_opportunity_fragment_ids
        )
        observed_fraction = observed_mutant_fragment_fraction(
            mutant_fragment_count,
            informative_fragment_count,
        )
        filter_reasons = _call_filter_reasons(
            mutant_fragment_count=mutant_fragment_count,
            informative_fragment_count=informative_fragment_count,
            observed_fraction=observed_fraction,
            partial_observation=representative.is_partial_observation,
            r1_mutant_count=r1_mutant_count,
            r2_mutant_count=r2_mutant_count,
            r1_opportunity_count=r1_opportunity_count,
            r2_opportunity_count=r2_opportunity_count,
            conflicting_fragment_count=len(consensus.conflicting_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
            wild_type_fragment_count=wild_type_fragment_count,
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            canonical_allele=key,
            mutant_fragment_count=mutant_fragment_count,
            wild_type_fragment_count=wild_type_fragment_count,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
            r1_mutant_count=r1_mutant_count,
            r2_mutant_count=r2_mutant_count,
            r1_opportunity_count=r1_opportunity_count,
            r2_opportunity_count=r2_opportunity_count,
            concordant_fragment_count=len(consensus.concordant_fragment_ids),
            single_mate_fragment_count=len(consensus.single_mate_fragment_ids),
            conflicting_fragment_count=len(consensus.conflicting_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
            not_informative_fragment_count=len(
                consensus.not_informative_fragment_ids
            ),
            consolidated_members=consolidated_members_by_key.get(key, ()),
        )
        calls.append(call)
        representatives.extend(
            _sorted_representatives(representative_map[key].values())
        )

    calls.sort(key=_sort_key)
    representatives.sort(key=_representative_sort_key)
    return calls, representatives


def call_exact_itds_with_representatives(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
    min_copied_segment_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
    consolidation: ITDConsolidationSettings = ITDConsolidationSettings(),
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call exact-match ITDs and retain one alignment per unique support pattern."""
    min_copied_segment_length = _resolved_min_copied_segment_length(
        min_insert_length,
        min_copied_segment_length,
    )
    alignments = list(alignments)
    (
        grouped_itds,
        representative_map,
        consensus_by_key,
        consolidated_members_by_key,
    ) = _collect_exact_itd_support(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        min_copied_segment_length=min_copied_segment_length,
        require_in_frame=require_in_frame,
        evidence_filter=evidence_filter,
        consolidation=consolidation,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for key, itds in grouped_itds.items():
        representative = _representative_itd(itds, key)
        consensus = consensus_by_key[key]
        mutant_fragment_count = len(consensus.mutant_fragment_ids)
        wild_type_fragment_count = len(consensus.wild_type_fragment_ids)
        informative_fragment_count = (
            mutant_fragment_count + wild_type_fragment_count
        )
        r1_mutant_count = len(consensus.r1_mutant_fragment_ids)
        r2_mutant_count = len(consensus.r2_mutant_fragment_ids)
        r1_opportunity_count = len(
            consensus.r1_opportunity_fragment_ids
        )
        r2_opportunity_count = len(
            consensus.r2_opportunity_fragment_ids
        )
        observed_fraction = observed_mutant_fragment_fraction(
            mutant_fragment_count,
            informative_fragment_count,
        )
        filter_reasons = _call_filter_reasons(
            mutant_fragment_count=mutant_fragment_count,
            informative_fragment_count=informative_fragment_count,
            observed_fraction=observed_fraction,
            partial_observation=representative.is_partial_observation,
            r1_mutant_count=r1_mutant_count,
            r2_mutant_count=r2_mutant_count,
            r1_opportunity_count=r1_opportunity_count,
            r2_opportunity_count=r2_opportunity_count,
            conflicting_fragment_count=len(consensus.conflicting_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
            wild_type_fragment_count=wild_type_fragment_count,
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            canonical_allele=key,
            mutant_fragment_count=mutant_fragment_count,
            wild_type_fragment_count=wild_type_fragment_count,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
            r1_mutant_count=r1_mutant_count,
            r2_mutant_count=r2_mutant_count,
            r1_opportunity_count=r1_opportunity_count,
            r2_opportunity_count=r2_opportunity_count,
            concordant_fragment_count=len(consensus.concordant_fragment_ids),
            single_mate_fragment_count=len(consensus.single_mate_fragment_ids),
            conflicting_fragment_count=len(consensus.conflicting_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
            not_informative_fragment_count=len(
                consensus.not_informative_fragment_ids
            ),
            consolidated_members=consolidated_members_by_key.get(key, ()),
        )
        calls.append(call)
        representatives.extend(
            _sorted_representatives(representative_map[key].values())
        )

    calls.sort(key=_sort_key)
    representatives.sort(key=_representative_sort_key)
    return calls, representatives


def _representative_itd(
    itds: list[ITD],
    canonical_allele: CanonicalInsertionAllele,
) -> ITD:
    return min(
        itds,
        key=lambda itd: (
            (
                itd.insertion.start != canonical_allele.start
                or itd.insertion.sequence != canonical_allele.sequence
            ),
            -itd.length,
            itd.spacer_length,
            itd.copied_segment_start,
            itd.insertion.start,
            itd.insertion.sequence,
            itd.insertion.read_id,
            itd.insertion.fragment_id,
            itd.insertion.direction,
        ),
    )


def _sort_key(call: ITDCall) -> tuple[int, int, str, str, str]:
    return (
        call.itd.insertion.start,
        call.itd.copied_segment_start,
        call.itd.spacer_prefix,
        call.itd.spacer_suffix,
        call.itd.insertion.sequence,
    )


def _support_signature(
    alignment: Alignment,
    itd: ITD,
    reference: str,
    *,
    flank_size: int,
) -> str:
    observed_bases = _observed_bases_by_reference_position(alignment)
    left_start = max(0, itd.insertion.start - flank_size + 1)
    left = "".join(
        observed_bases.get(position, "-")
        for position in range(left_start, itd.insertion.start + 1)
    )
    right_end = min(len(reference), itd.copied_segment_start + flank_size)
    right = "".join(
        observed_bases.get(position, "-")
        for position in range(itd.copied_segment_start, right_end)
    )
    return f"{left}[{itd.copied_segment_sequence}]{right}"


def _observed_bases_by_reference_position(
    alignment: Alignment,
) -> dict[int, str]:
    ref_pos = -1
    observed_bases: dict[int, str] = {}

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-":
                observed_bases[ref_pos] = read_base

    return observed_bases


def _collect_exact_itd_support(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int,
    min_copied_segment_length: int,
    require_in_frame: bool,
    evidence_filter: InsertionEvidenceFilter | None,
    consolidation: ITDConsolidationSettings,
) -> tuple[
    dict[ITDCallKey, list[ITD]],
    SupportRepresentativeMap,
    dict[ITDCallKey, FragmentConsensusSupport],
    dict[ITDCallKey, tuple[ConsolidatedAlleleMember, ...]],
]:
    alignments = list(alignments)
    grouped_itds: dict[ITDCallKey, list[ITD]] = defaultdict(list)
    representative_map: dict[
        ITDCallKey,
        dict[str, tuple[ITD, Alignment, bool]],
    ] = defaultdict(dict)
    fragment_ids_by_signature: dict[ITDCallKey, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    insert_sequences_by_fragment: dict[
        ITDCallKey, dict[str, dict[str, int]]
    ] = defaultdict(lambda: defaultdict(dict))

    observations: list[_CandidateObservation] = []
    for alignment in alignments:
        insertions = extract_insertions(
            alignment,
            min_length=min_insert_length,
            require_in_frame=require_in_frame,
            evidence_filter=None,
        )
        passing_insertions = set(
            extract_insertions(
                alignment,
                min_length=min_insert_length,
                require_in_frame=require_in_frame,
                evidence_filter=evidence_filter,
            )
        )
        for insertion in insertions:
            itd = classify_exact_itd(
                insertion,
                reference,
                min_copied_segment_length=min_copied_segment_length,
            )
            if itd is None:
                continue
            observations.append(
                _CandidateObservation(
                    itd=itd,
                    alignment=alignment,
                    observed_allele=canonicalize_insertion_allele(
                        itd.insertion,
                        reference,
                    ),
                    canonical_allele=canonicalize_insertion_allele(
                        itd.insertion,
                        reference,
                    ),
                    passes_evidence=insertion in passing_insertions,
                )
            )

    observations, consolidated_members_by_key = _consolidate_observations(
        observations,
        reference,
        consolidation,
    )
    consensus_by_key = _fragment_consensus_support(
        observations,
        alignments,
        reference=reference,
        evidence_filter=evidence_filter,
    )
    for observation in observations:
        itd = observation.itd
        alignment = observation.alignment

        key = observation.canonical_allele
        expected_sequence = _expected_insertion_sequence(itd)
        grouped_itds[key].append(itd)
        signature = _support_signature(
            alignment,
            itd,
            reference,
            flank_size=itd.length,
        )
        _set_best_representative(
            representative_map[key],
            signature,
            itd,
            alignment,
            observation.passes_evidence,
        )
        if (
            observation.passes_evidence
            and alignment.fragment_id
            in consensus_by_key[key].mutant_fragment_ids
        ):
            fragment_ids_by_signature[key][signature].add(alignment.fragment_id)
            insert_sequences_by_fragment[key][alignment.fragment_id][
                itd.insertion.sequence
            ] = _sequence_mismatches(itd.insertion.sequence, expected_sequence)

    finalized_map: SupportRepresentativeMap = defaultdict(dict)
    for key, alignments_by_signature in representative_map.items():
        insert_sequence_supports = _insert_sequence_supports(
            insert_sequences_by_fragment[key]
        )
        for signature, (
            itd,
            alignment,
            _,
        ) in alignments_by_signature.items():
            expected_sequence = _expected_insertion_sequence(itd)
            finalized_map[key][signature] = UniqueSupportRepresentative(
                itd=itd,
                signature=signature,
                alignment=alignment,
                support_count=len(fragment_ids_by_signature[key][signature]),
                exact_support_count=len(fragment_ids_by_signature[key][signature]),
                mismatches=_sequence_mismatches(
                    itd.insertion.sequence,
                    expected_sequence,
                ),
                insert_sequence_supports=insert_sequence_supports,
                canonical_allele=key,
            )

    return (
        grouped_itds,
        finalized_map,
        consensus_by_key,
        consolidated_members_by_key,
    )


def _collect_fuzzy_itd_support(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_copy_mismatch_rate: float,
    min_insert_length: int,
    min_copied_segment_length: int,
    require_in_frame: bool,
    evidence_filter: InsertionEvidenceFilter | None,
    consolidation: ITDConsolidationSettings,
) -> tuple[
    dict[ITDCallKey, list[ITD]],
    SupportRepresentativeMap,
    dict[ITDCallKey, FragmentConsensusSupport],
    dict[ITDCallKey, tuple[ConsolidatedAlleleMember, ...]],
]:
    alignments = list(alignments)
    grouped_itds: dict[ITDCallKey, list[ITD]] = defaultdict(list)
    representative_map: dict[
        ITDCallKey,
        dict[str, tuple[ITD, Alignment, bool]],
    ] = defaultdict(dict)
    fragment_ids_by_signature: dict[ITDCallKey, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    exact_fragment_ids_by_signature: dict[
        ITDCallKey, dict[str, set[str]]
    ] = defaultdict(lambda: defaultdict(set))
    fuzzy_example_sequence_by_signature: dict[
        ITDCallKey, dict[str, str]
    ] = defaultdict(dict)
    insert_sequences_by_fragment: dict[
        ITDCallKey, dict[str, dict[str, int]]
    ] = defaultdict(lambda: defaultdict(dict))

    observations: list[_CandidateObservation] = []
    for alignment in alignments:
        insertions = extract_insertions(
            alignment,
            min_length=min_insert_length,
            require_in_frame=require_in_frame,
            evidence_filter=None,
        )
        passing_insertions = set(
            extract_insertions(
                alignment,
                min_length=min_insert_length,
                require_in_frame=require_in_frame,
                evidence_filter=evidence_filter,
            )
        )
        for insertion in insertions:
            itd = classify_fuzzy_itd(
                insertion,
                reference,
                max_copy_mismatch_rate=max_copy_mismatch_rate,
                min_copied_segment_length=min_copied_segment_length,
            )
            if itd is None:
                continue

            observations.append(
                _CandidateObservation(
                    itd=itd,
                    alignment=alignment,
                    observed_allele=canonicalize_insertion_allele(
                        insertion,
                        reference,
                    ),
                    canonical_allele=canonicalize_insertion(
                        start=insertion.start,
                        sequence=_expected_insertion_sequence(itd),
                        reference=reference,
                        trailing=insertion.trailing,
                    ),
                    passes_evidence=insertion in passing_insertions,
                )
            )

    observations, consolidated_members_by_key = _consolidate_observations(
        observations,
        reference,
        consolidation,
    )
    consensus_by_key = _fragment_consensus_support(
        observations,
        alignments,
        reference=reference,
        evidence_filter=evidence_filter,
    )
    for observation in observations:
        itd = observation.itd
        alignment = observation.alignment
        insertion = itd.insertion

        key = observation.canonical_allele
        expected_sequence = _expected_insertion_sequence(itd)
        sequence_mismatches = _sequence_mismatches(
            insertion.sequence,
            expected_sequence,
        )
        grouped_itds[key].append(itd)
        signature = _support_signature(
            alignment,
            itd,
            reference,
            flank_size=itd.length,
        )
        _set_best_representative(
            representative_map[key],
            signature,
            itd,
            alignment,
            observation.passes_evidence,
        )
        is_supporting = (
            observation.passes_evidence
            and alignment.fragment_id
            in consensus_by_key[key].mutant_fragment_ids
        )
        if is_supporting:
            fragment_ids_by_signature[key][signature].add(alignment.fragment_id)
            insert_sequences_by_fragment[key][alignment.fragment_id][
                insertion.sequence
            ] = sequence_mismatches

        if is_supporting and sequence_mismatches == 0:
            exact_fragment_ids_by_signature[key][signature].add(
                alignment.fragment_id
            )
            continue

        if is_supporting:
            fuzzy_example_sequence_by_signature[key].setdefault(
                signature,
                insertion.sequence,
            )

    finalized_map: SupportRepresentativeMap = defaultdict(dict)
    for key, alignments_by_signature in representative_map.items():
        insert_sequence_supports = _insert_sequence_supports(
            insert_sequences_by_fragment[key]
        )
        for signature, (
            itd,
            alignment,
            _,
        ) in alignments_by_signature.items():
            fragment_ids = fragment_ids_by_signature[key][signature]
            exact_fragment_ids = exact_fragment_ids_by_signature[key][signature]
            fuzzy_only_count = len(fragment_ids - exact_fragment_ids)
            expected_sequence = _expected_insertion_sequence(itd)
            finalized_map[key][signature] = UniqueSupportRepresentative(
                itd=itd,
                signature=signature,
                alignment=alignment,
                support_count=len(fragment_ids),
                exact_support_count=len(exact_fragment_ids),
                fuzzy_only_support_count=fuzzy_only_count,
                fuzzy_example_sequence=fuzzy_example_sequence_by_signature[key].get(
                    signature
                ),
                mismatches=_sequence_mismatches(
                    itd.insertion.sequence,
                    expected_sequence,
                ),
                insert_sequence_supports=insert_sequence_supports,
                canonical_allele=key,
            )

    return (
        grouped_itds,
        finalized_map,
        consensus_by_key,
        consolidated_members_by_key,
    )


def _consolidate_observations(
    observations: list[_CandidateObservation],
    reference: str,
    settings: ITDConsolidationSettings,
) -> tuple[
    list[_CandidateObservation],
    dict[ITDCallKey, tuple[ConsolidatedAlleleMember, ...]],
]:
    """Assign weak compatible alleles directly to dominant anchor alleles."""
    if not settings.enabled:
        return observations, {}

    fragment_ids_by_key: dict[ITDCallKey, set[str]] = defaultdict(set)
    keys: set[ITDCallKey] = set()
    for observation in observations:
        key = observation.canonical_allele
        keys.add(key)
        if observation.passes_evidence:
            fragment_ids_by_key[key].add(observation.alignment.fragment_id)

    support_by_key = {
        key: len(fragment_ids_by_key.get(key, set()))
        for key in keys
    }
    ordered_keys = sorted(
        keys,
        key=lambda key: (
            -support_by_key[key],
            key.start,
            key.sequence,
            key.trailing,
        ),
    )
    anchors: list[ITDCallKey] = []
    assignment: dict[ITDCallKey, ITDCallKey] = {}
    member_details: dict[
        ITDCallKey,
        list[ConsolidatedAlleleMember],
    ] = defaultdict(list)

    for key in ordered_keys:
        compatible = [
            (
                _consolidation_rank(
                    key,
                    anchor,
                    support_by_key[anchor],
                    reference,
                ),
                anchor,
            )
            for anchor in anchors
            if _can_consolidate_allele(
                key,
                anchor,
                minor_support=support_by_key[key],
                anchor_support=support_by_key[anchor],
                reference=reference,
                settings=settings,
            )
        ]
        compatible.sort(key=lambda item: (item[0], item[1]))
        if compatible and (
            len(compatible) == 1
            or compatible[0][0] != compatible[1][0]
        ):
            rank, anchor = compatible[0]
            assignment[key] = anchor
            (
                allele_mismatch_rate,
                breakpoint_shift_rate,
                _,
            ) = rank
            allele_mismatches = _alternate_sequence_mismatches(
                key,
                anchor,
                reference,
            )
            breakpoint_shift = abs(key.start - anchor.start)
            member_details[anchor].append(
                ConsolidatedAlleleMember(
                    allele=key,
                    fragment_count=support_by_key[key],
                    allele_mismatches=allele_mismatches,
                    allele_mismatch_rate=allele_mismatch_rate,
                    breakpoint_shift=breakpoint_shift,
                    breakpoint_shift_rate=breakpoint_shift_rate,
                    reason=(
                        "same-breakpoint sequence error"
                        if breakpoint_shift == 0
                        else "nearby-breakpoint local-haplotype match"
                    ),
                )
            )
            continue

        assignment[key] = key
        if support_by_key[key] >= settings.min_anchor_fragment_count:
            anchors.append(key)

    consolidated = [
        replace(
            observation,
            canonical_allele=assignment[observation.canonical_allele],
            consolidated_from=(
                observation.canonical_allele
                if assignment[observation.canonical_allele]
                != observation.canonical_allele
                else None
            ),
        )
        for observation in observations
    ]
    finalized_details = {
        anchor: tuple(
            sorted(
                members,
                key=lambda member: (
                    member.allele_mismatch_rate,
                    member.breakpoint_shift_rate,
                    -member.fragment_count,
                    member.allele,
                ),
            )
        )
        for anchor, members in member_details.items()
    }
    return consolidated, finalized_details


def _can_consolidate_allele(
    minor: ITDCallKey,
    anchor: ITDCallKey,
    *,
    minor_support: int,
    anchor_support: int,
    reference: str,
    settings: ITDConsolidationSettings,
) -> bool:
    if minor_support < 1:
        return False
    if minor.trailing or anchor.trailing:
        return False
    if len(minor.sequence) != len(anchor.sequence):
        return False
    allele_length = len(minor.sequence)
    breakpoint_shift = abs(minor.start - anchor.start)
    if (
        breakpoint_shift / allele_length
        > settings.max_breakpoint_shift_rate
    ):
        return False
    if anchor_support < settings.min_anchor_fragment_count:
        return False
    if minor_support > (
        anchor_support * settings.max_minor_to_anchor_support_ratio
    ):
        return False
    return (
        _alternate_sequence_mismatches(minor, anchor, reference)
        / allele_length
        <= settings.max_allele_mismatch_rate
    )


def _consolidation_rank(
    minor: ITDCallKey,
    anchor: ITDCallKey,
    anchor_support: int,
    reference: str,
) -> tuple[float, float, int]:
    allele_length = len(minor.sequence)
    return (
        _alternate_sequence_mismatches(minor, anchor, reference)
        / allele_length,
        abs(minor.start - anchor.start) / allele_length,
        -anchor_support,
    )


def _alternate_sequence_mismatches(
    first: ITDCallKey,
    second: ITDCallKey,
    reference: str,
) -> int:
    first_alt = first.alternate_sequence(reference)
    second_alt = second.alternate_sequence(reference)
    if len(first_alt) != len(second_alt):
        return max(len(first_alt), len(second_alt))
    return _sequence_mismatches(first_alt, second_alt)


def _fragment_consensus_support(
    observations: Iterable[_CandidateObservation],
    alignments: Iterable[Alignment],
    *,
    reference: str,
    evidence_filter: InsertionEvidenceFilter | None,
) -> dict[ITDCallKey, FragmentConsensusSupport]:
    """Assign every fragment one candidate-specific evidence state."""
    observations = list(observations)
    alignments = list(alignments)
    candidates_by_fragment_direction: dict[
        str,
        dict[
            str,
            set[
                tuple[
                    ITDCallKey,
                    CanonicalInsertionAllele,
                    CanonicalInsertionAllele | None,
                ]
            ],
        ],
    ] = defaultdict(lambda: defaultdict(set))
    rejected_candidate_keys_by_fragment: dict[str, set[ITDCallKey]] = defaultdict(
        set
    )
    all_keys: set[ITDCallKey] = set()
    for observation in observations:
        key = observation.canonical_allele
        all_keys.add(key)
        if not observation.passes_evidence:
            rejected_candidate_keys_by_fragment[
                observation.alignment.fragment_id
            ].add(key)
            continue
        candidates_by_fragment_direction[observation.alignment.fragment_id][
            observation.alignment.direction
        ].add(
            (
                key,
                observation.observed_allele,
                observation.consolidated_from,
            )
        )

    alignments_by_fragment_direction: dict[
        str, dict[str, list[Alignment]]
    ] = defaultdict(lambda: defaultdict(list))
    for alignment in alignments:
        alignments_by_fragment_direction[alignment.fragment_id][
            alignment.direction
        ].append(alignment)

    category_sets: dict[ITDCallKey, dict[str, set[str]]] = {
        key: {
            "supporting": set(),
            "wild_type": set(),
            "concordant": set(),
            "single_mate": set(),
            "conflicting": set(),
            "unresolved": set(),
            "not_informative": set(),
            "r1_mutant": set(),
            "r2_mutant": set(),
            "r1_opportunity": set(),
            "r2_opportunity": set(),
        }
        for key in all_keys
    }

    for (
        fragment_id,
        candidates_by_direction,
    ) in candidates_by_fragment_direction.items():
        forward_candidates = candidates_by_direction.get("forward", set())
        reverse_candidates = candidates_by_direction.get("reverse", set())
        forward_keys = {key for key, _, _ in forward_candidates}
        reverse_keys = {key for key, _, _ in reverse_candidates}
        candidate_keys = forward_keys | reverse_keys
        if len(forward_keys) > 1 or len(reverse_keys) > 1:
            for key in candidate_keys:
                category_sets[key]["unresolved"].add(fragment_id)
            continue
        forward_observed = {
            observed for _, observed, _ in forward_candidates
        }
        reverse_observed = {
            observed for _, observed, _ in reverse_candidates
        }
        forward_unconsolidated_observed = {
            observed
            for _, observed, consolidated_from in forward_candidates
            if consolidated_from is None
        }
        reverse_unconsolidated_observed = {
            observed
            for _, observed, consolidated_from in reverse_candidates
            if consolidated_from is None
        }
        if (
            len(forward_observed) > 1
            and len(forward_unconsolidated_observed) > 1
        ) or (
            len(reverse_observed) > 1
            and len(reverse_unconsolidated_observed) > 1
        ):
            for key in candidate_keys:
                category_sets[key]["unresolved"].add(fragment_id)
            continue

        if forward_candidates and reverse_candidates:
            forward_key = next(iter(forward_keys))
            reverse_key = next(iter(reverse_keys))
            if forward_key != reverse_key:
                for key in candidate_keys:
                    category_sets[key]["conflicting"].add(fragment_id)
                continue
            if (
                forward_observed != reverse_observed
                and len(
                    forward_unconsolidated_observed
                    | reverse_unconsolidated_observed
                )
                > 1
            ):
                category_sets[forward_key]["unresolved"].add(fragment_id)
                continue
            key = forward_key
            category_sets[key]["supporting"].add(fragment_id)
            category_sets[key]["concordant"].add(fragment_id)
            continue

        key = next(iter(forward_keys or reverse_keys))
        candidate_direction = "forward" if forward_candidates else "reverse"
        opposite_direction = (
            "reverse" if candidate_direction == "forward" else "forward"
        )
        opposite_supports_wild_type = any(
            _alignment_supports_wild_type_junction(
                alignment,
                key.start,
                reference,
                evidence_filter,
            )
            for alignment in alignments_by_fragment_direction[fragment_id].get(
                opposite_direction,
                [],
            )
        )
        if opposite_supports_wild_type:
            category_sets[key]["conflicting"].add(fragment_id)
        else:
            category_sets[key]["supporting"].add(fragment_id)
            category_sets[key]["single_mate"].add(fragment_id)

    all_fragment_ids = set(alignments_by_fragment_direction)
    evidence_states = (
        "supporting",
        "wild_type",
        "conflicting",
        "unresolved",
        "not_informative",
    )
    for key in all_keys:
        for fragment_id in all_fragment_ids:
            categories = category_sets[key]
            if any(fragment_id in categories[state] for state in evidence_states):
                continue
            fragment_alignments = [
                alignment
                for direction_alignments in alignments_by_fragment_direction[
                    fragment_id
                ].values()
                for alignment in direction_alignments
            ]
            if key in rejected_candidate_keys_by_fragment.get(fragment_id, set()):
                categories["unresolved"].add(fragment_id)
            elif any(
                _alignment_supports_wild_type_junction(
                    alignment,
                    key.start,
                    reference,
                    evidence_filter,
                )
                for alignment in fragment_alignments
            ):
                categories["wild_type"].add(fragment_id)
            elif any(
                spans_insertion_site(alignment, key.start)
                for alignment in fragment_alignments
            ):
                categories["unresolved"].add(fragment_id)
            else:
                categories["not_informative"].add(fragment_id)

    for key, categories in category_sets.items():
        informative_fragment_ids = (
            categories["supporting"] | categories["wild_type"]
        )
        for fragment_id in informative_fragment_ids:
            for direction in ("forward", "reverse"):
                read = "r1" if direction == "forward" else "r2"
                direction_candidates = candidates_by_fragment_direction[
                    fragment_id
                ].get(direction, set())
                supports_candidate = any(
                    candidate_key == key
                    for candidate_key, _, _ in direction_candidates
                )
                direction_alignments = alignments_by_fragment_direction[
                    fragment_id
                ].get(direction, [])
                supports_wild_type = any(
                    _alignment_supports_wild_type_junction(
                        alignment,
                        key.start,
                        reference,
                        evidence_filter,
                    )
                    for alignment in direction_alignments
                )
                if supports_candidate:
                    categories[f"{read}_mutant"].add(fragment_id)
                if supports_candidate or supports_wild_type:
                    categories[f"{read}_opportunity"].add(fragment_id)

    return {
        key: FragmentConsensusSupport(
            mutant_fragment_ids=frozenset(categories["supporting"]),
            wild_type_fragment_ids=frozenset(categories["wild_type"]),
            concordant_fragment_ids=frozenset(categories["concordant"]),
            single_mate_fragment_ids=frozenset(categories["single_mate"]),
            conflicting_fragment_ids=frozenset(categories["conflicting"]),
            unresolved_fragment_ids=frozenset(categories["unresolved"]),
            not_informative_fragment_ids=frozenset(
                categories["not_informative"]
            ),
            r1_mutant_fragment_ids=frozenset(categories["r1_mutant"]),
            r2_mutant_fragment_ids=frozenset(categories["r2_mutant"]),
            r1_opportunity_fragment_ids=frozenset(
                categories["r1_opportunity"]
            ),
            r2_opportunity_fragment_ids=frozenset(
                categories["r2_opportunity"]
            ),
        )
        for key, categories in category_sets.items()
    }


def _alignment_supports_wild_type_junction(
    alignment: Alignment,
    site: int,
    reference: str,
    evidence_filter: InsertionEvidenceFilter | None,
) -> bool:
    """Return whether an alignment gives high-quality WT junction evidence."""
    if not spans_insertion_site(alignment, site):
        return False

    flank = evidence_filter.junction_flank_size if evidence_filter else 1
    left_start = site - flank + 1
    right_end = site + flank
    if left_start < 0 or right_end >= len(reference):
        return False
    required_positions = list(range(left_start, site + 1)) + list(
        range(site + 1, right_end + 1)
    )

    columns_by_position: dict[int, int] = {}
    ref_pos = -1
    for column, ref_base in enumerate(alignment.aligned_reference):
        if ref_base != "-":
            ref_pos += 1
            columns_by_position[ref_pos] = column

    try:
        columns = [columns_by_position[position] for position in required_positions]
    except KeyError:
        return False
    if any(right != left + 1 for left, right in zip(columns, columns[1:])):
        return False

    for position, column in zip(required_positions, columns, strict=True):
        read_base = alignment.aligned_read[column]
        if read_base != reference[position]:
            return False
        if evidence_filter is not None and alignment.aligned_qualities:
            quality = alignment.aligned_qualities[column]
            if (
                quality is None
                or quality < evidence_filter.min_junction_anchor_quality
            ):
                return False
    return True


def _set_best_representative(
    representatives: dict[str, tuple[ITD, Alignment, bool]],
    signature: str,
    itd: ITD,
    alignment: Alignment,
    passes_evidence: bool,
) -> None:
    current = representatives.get(signature)
    candidate_key = (
        not passes_evidence,
        _sequence_mismatches(
            itd.insertion.sequence,
            _expected_insertion_sequence(itd),
        ),
        alignment.read_id,
    )
    if current is None:
        representatives[signature] = (itd, alignment, passes_evidence)
        return

    current_itd, current_alignment, current_passes_evidence = current
    current_key = (
        not current_passes_evidence,
        _sequence_mismatches(
            current_itd.insertion.sequence,
            _expected_insertion_sequence(current_itd),
        ),
        current_alignment.read_id,
    )
    if candidate_key < current_key:
        representatives[signature] = (itd, alignment, passes_evidence)


def _insert_sequence_supports(
    sequences_by_fragment: dict[str, dict[str, int]],
) -> tuple[InsertSequenceSupport, ...]:
    fragment_ids_by_sequence: dict[str, set[str]] = defaultdict(set)
    mismatches_by_sequence: dict[str, int] = {}

    for fragment_id, sequence_mismatches in sequences_by_fragment.items():
        sequence, mismatches = min(
            sequence_mismatches.items(),
            key=lambda item: (item[1], item[0]),
        )
        fragment_ids_by_sequence[sequence].add(fragment_id)
        mismatches_by_sequence[sequence] = mismatches

    supports = [
        InsertSequenceSupport(
            sequence=sequence,
            support_count=len(fragment_ids),
            mismatches=mismatches_by_sequence[sequence],
        )
        for sequence, fragment_ids in fragment_ids_by_sequence.items()
    ]
    return tuple(
        sorted(
            supports,
            key=lambda support: (
                support.mismatches,
                -support.support_count,
                support.sequence,
            ),
        )
    )


def _sequence_mismatches(observed: str, expected: str) -> int:
    return sum(
        1
        for observed_base, expected_base in zip(observed, expected, strict=True)
        if observed_base != expected_base
    )


def _sorted_representatives(
    representatives: Iterable[UniqueSupportRepresentative],
) -> list[UniqueSupportRepresentative]:
    return sorted(representatives, key=_representative_sort_key)


def _representative_sort_key(
    representative: UniqueSupportRepresentative,
) -> tuple[int, int, str, str, str, int, str, str]:
    return (
        representative.itd.insertion.start,
        representative.itd.copied_segment_start,
        representative.itd.copied_segment_sequence,
        representative.itd.spacer_prefix,
        representative.itd.spacer_suffix,
        -representative.support_count,
        representative.signature,
        representative.alignment.read_id,
    )


def _call_filter_reasons(
    *,
    mutant_fragment_count: int,
    informative_fragment_count: int,
    observed_fraction: float,
    partial_observation: bool,
    r1_mutant_count: int,
    r2_mutant_count: int,
    r1_opportunity_count: int,
    r2_opportunity_count: int,
    conflicting_fragment_count: int,
    unresolved_fragment_count: int,
    wild_type_fragment_count: int,
    filters: ITDFilter,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if partial_observation:
        reasons.append("PARTIAL_OBSERVATION")
    if mutant_fragment_count == 0 and conflicting_fragment_count:
        reasons.append("ONLY_CONFLICTING_MATE_EVIDENCE")
    if mutant_fragment_count == 0 and unresolved_fragment_count:
        reasons.append("ONLY_UNRESOLVED_EVIDENCE")
    if (
        conflicting_fragment_count + unresolved_fragment_count
        > mutant_fragment_count + wild_type_fragment_count
    ):
        reasons.append("AMBIGUOUS_EVIDENCE_DOMINATES")
    if mutant_fragment_count < filters.min_mutant_fragment_count:
        reasons.append("LOW_MUTANT_FRAGMENT_COUNT")
    if informative_fragment_count < filters.min_informative_fragment_count:
        reasons.append("LOW_INFORMATIVE_FRAGMENT_COUNT")
    if observed_fraction < filters.min_observed_mutant_fragment_fraction:
        reasons.append("LOW_MUTANT_FRAGMENT_FRACTION")
    if (
        r1_opportunity_count >= filters.min_directional_opportunities
        and r2_opportunity_count >= filters.min_directional_opportunities
    ):
        r1_mutant_fraction = (
            r1_mutant_count / r1_opportunity_count
        )
        r2_mutant_fraction = (
            r2_mutant_count / r2_opportunity_count
        )
        summed_directional_fractions = (
            r1_mutant_fraction + r2_mutant_fraction
        )
        if summed_directional_fractions > 0 and (
            max(r1_mutant_fraction, r2_mutant_fraction)
            / summed_directional_fractions
            > filters.max_directional_mutant_fraction_share
        ):
            reasons.append("DIRECTION_BIAS")
    return tuple(reasons)


def _expected_insertion_sequence(itd: ITD) -> str:
    return f"{itd.spacer_prefix}{itd.copied_segment_sequence}{itd.spacer_suffix}"


def _resolved_min_copied_segment_length(
    min_insert_length: int,
    min_copied_segment_length: int | None,
) -> int:
    if min_insert_length < 1:
        raise ValueError("min_insert_length must be at least 1")
    resolved = (
        min_insert_length
        if min_copied_segment_length is None
        else min_copied_segment_length
    )
    if resolved < 1:
        raise ValueError("min_copied_segment_length must be at least 1")
    return resolved
