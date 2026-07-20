"""Build reportable ITD calls from aligned reads."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
import math
import warnings

from .coverage import (
    interbase_fragment_ids,
    observed_supporting_fragment_fraction,
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
    """An ITD with transparent post-filter fragment counts and fraction."""

    itd: ITD
    supporting_fragment_count: int
    spanning_fragment_count: int
    observed_supporting_fragment_fraction: float
    status: str = "PASS"
    filter_reasons: tuple[str, ...] = ()
    forward_support_count: int = field(default=0, compare=False)
    reverse_support_count: int = field(default=0, compare=False)
    concordant_fragment_count: int = field(default=0, compare=False)
    single_mate_fragment_count: int = field(default=0, compare=False)
    discordant_fragment_count: int = field(default=0, compare=False)
    unresolved_fragment_count: int = field(default=0, compare=False)

    def __post_init__(self) -> None:
        expected_fraction = observed_supporting_fragment_fraction(
            self.supporting_fragment_count,
            self.spanning_fragment_count,
        )
        if not math.isclose(
            self.observed_supporting_fragment_fraction,
            expected_fraction,
            rel_tol=0,
            abs_tol=1e-12,
        ):
            raise ValueError(
                "observed_supporting_fragment_fraction must equal "
                "supporting_fragment_count / spanning_fragment_count"
            )

    @property
    def passes_filters(self) -> bool:
        """Return whether the call passes the configured thresholds."""
        return self.status == "PASS"

    @property
    def support_count(self) -> int:
        """Deprecated alias for ``supporting_fragment_count``."""
        _warn_deprecated_call_attribute("support_count", "supporting_fragment_count")
        return self.supporting_fragment_count

    @property
    def coverage(self) -> int:
        """Deprecated alias for ``spanning_fragment_count``."""
        _warn_deprecated_call_attribute("coverage", "spanning_fragment_count")
        return self.spanning_fragment_count

    @property
    def vaf(self) -> float:
        """Deprecated alias; this quantity is not a validated VAF."""
        _warn_deprecated_call_attribute(
            "vaf",
            "observed_supporting_fragment_fraction",
        )
        return self.observed_supporting_fragment_fraction


@dataclass(frozen=True)
class ITDFilter:
    """Thresholds used to label exact-match ITD calls."""

    min_supporting_fragment_count: int = 3
    min_spanning_fragment_count: int = 10
    min_observed_supporting_fragment_fraction: float = 0.01
    max_single_direction_fraction: float = 0.90
    min_directional_observations: int = 5

    def __post_init__(self) -> None:
        if self.min_supporting_fragment_count < 1:
            raise ValueError("min_supporting_fragment_count must be at least 1")
        if self.min_spanning_fragment_count < 0:
            raise ValueError("min_spanning_fragment_count must not be negative")
        if not 0 <= self.min_observed_supporting_fragment_fraction <= 1:
            raise ValueError(
                "min_observed_supporting_fragment_fraction must be between 0 and 1"
            )
        if not 0.5 <= self.max_single_direction_fraction <= 1:
            raise ValueError(
                "max_single_direction_fraction must be between 0.5 and 1"
            )
        if self.min_directional_observations < 1:
            raise ValueError("min_directional_observations must be at least 1")

    @property
    def min_support_count(self) -> int:
        """Deprecated alias for ``min_supporting_fragment_count``."""
        _warn_deprecated_call_attribute(
            "min_support_count",
            "min_supporting_fragment_count",
        )
        return self.min_supporting_fragment_count

    @property
    def min_coverage(self) -> int:
        """Deprecated alias for ``min_spanning_fragment_count``."""
        _warn_deprecated_call_attribute(
            "min_coverage",
            "min_spanning_fragment_count",
        )
        return self.min_spanning_fragment_count

    @property
    def min_vaf(self) -> float:
        """Deprecated alias for the minimum observed fragment fraction."""
        _warn_deprecated_call_attribute(
            "min_vaf",
            "min_observed_supporting_fragment_fraction",
        )
        return self.min_observed_supporting_fragment_fraction


@dataclass(frozen=True)
class UniqueSupportRepresentative:
    """One representative alignment for a unique local ITD support pattern."""

    itd: ITD
    signature: str
    alignment: Alignment
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
    """Fragment eligibility and mate-consensus categories for one call."""

    supporting_fragment_ids: frozenset[str] = frozenset()
    concordant_fragment_ids: frozenset[str] = frozenset()
    single_mate_fragment_ids: frozenset[str] = frozenset()
    discordant_fragment_ids: frozenset[str] = frozenset()
    unresolved_fragment_ids: frozenset[str] = frozenset()


@dataclass(frozen=True)
class _CandidateObservation:
    itd: ITD
    alignment: Alignment


# The breakpoint is part of the event identity.  Calls with the same copied
# tract but different insertion sites need separate support and spanning counts.
ITDCallKey = tuple[int, int, str, str, str, bool]
SupportRepresentativeMap = dict[ITDCallKey, dict[str, UniqueSupportRepresentative]]


def call_exact_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
    min_tandem_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
) -> list[ITDCall]:
    """Call exact-match ITDs and attach fragment counts and observed fraction."""
    calls, _ = call_exact_itds_with_representatives(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        min_tandem_length=min_tandem_length,
        require_in_frame=require_in_frame,
        filters=filters,
        evidence_filter=evidence_filter,
    )
    return calls


def call_fuzzy_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_mismatches: int,
    min_insert_length: int = 6,
    min_tandem_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
) -> list[ITDCall]:
    """Call fuzzy-match ITDs and attach fragment counts and observed fraction."""
    if max_mismatches < 0:
        raise ValueError("max_mismatches must not be negative")
    calls, _ = call_fuzzy_itds_with_representatives(
        alignments,
        reference,
        max_mismatches=max_mismatches,
        min_insert_length=min_insert_length,
        min_tandem_length=min_tandem_length,
        require_in_frame=require_in_frame,
        filters=filters,
        evidence_filter=evidence_filter,
    )
    return calls


def call_fuzzy_itds_with_representatives(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_mismatches: int,
    min_insert_length: int = 6,
    min_tandem_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call fuzzy-match ITDs and retain one alignment per unique support pattern."""
    if max_mismatches < 0:
        raise ValueError("max_mismatches must not be negative")
    min_tandem_length = _resolved_min_tandem_length(
        min_insert_length,
        min_tandem_length,
    )
    alignments = list(alignments)
    spanning_fragments_by_site = interbase_fragment_ids(alignments)
    grouped_itds, representative_map, consensus_by_key = _collect_fuzzy_itd_support(
        alignments,
        reference,
        max_mismatches=max_mismatches,
        min_insert_length=min_insert_length,
        min_tandem_length=min_tandem_length,
        require_in_frame=require_in_frame,
        evidence_filter=evidence_filter,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for itds in grouped_itds.values():
        representative = _representative_itd(itds)
        key = _itd_call_key(representative)
        consensus = consensus_by_key[key]
        supporting_fragment_count = len(consensus.supporting_fragment_ids)
        forward_support_count, reverse_support_count = _direction_support_counts(
            itds,
            consensus.supporting_fragment_ids,
        )
        eligible_spanning_fragments = spanning_fragments_by_site.get(
            representative.insertion.start,
            frozenset(),
        ) - consensus.discordant_fragment_ids - consensus.unresolved_fragment_ids
        spanning_fragment_count = len(eligible_spanning_fragments)
        observed_fraction = observed_supporting_fragment_fraction(
            supporting_fragment_count,
            spanning_fragment_count,
        )
        filter_reasons = _call_filter_reasons(
            supporting_fragment_count=supporting_fragment_count,
            spanning_fragment_count=spanning_fragment_count,
            observed_fraction=observed_fraction,
            partial_observation=representative.is_partial_observation,
            forward_support_count=forward_support_count,
            reverse_support_count=reverse_support_count,
            discordant_fragment_count=len(consensus.discordant_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            supporting_fragment_count=supporting_fragment_count,
            spanning_fragment_count=spanning_fragment_count,
            observed_supporting_fragment_fraction=observed_fraction,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
            forward_support_count=forward_support_count,
            reverse_support_count=reverse_support_count,
            concordant_fragment_count=len(consensus.concordant_fragment_ids),
            single_mate_fragment_count=len(consensus.single_mate_fragment_ids),
            discordant_fragment_count=len(consensus.discordant_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
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
    min_tandem_length: int | None = None,
    require_in_frame: bool = True,
    filters: ITDFilter = ITDFilter(),
    evidence_filter: InsertionEvidenceFilter | None = None,
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call exact-match ITDs and retain one alignment per unique support pattern."""
    min_tandem_length = _resolved_min_tandem_length(
        min_insert_length,
        min_tandem_length,
    )
    alignments = list(alignments)
    spanning_fragments_by_site = interbase_fragment_ids(alignments)
    grouped_itds, representative_map, consensus_by_key = _collect_exact_itd_support(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        min_tandem_length=min_tandem_length,
        require_in_frame=require_in_frame,
        evidence_filter=evidence_filter,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for itds in grouped_itds.values():
        representative = _representative_itd(itds)
        key = _itd_call_key(representative)
        consensus = consensus_by_key[key]
        supporting_fragment_count = len(consensus.supporting_fragment_ids)
        forward_support_count, reverse_support_count = _direction_support_counts(
            itds,
            consensus.supporting_fragment_ids,
        )
        eligible_spanning_fragments = spanning_fragments_by_site.get(
            representative.insertion.start,
            frozenset(),
        ) - consensus.discordant_fragment_ids - consensus.unresolved_fragment_ids
        spanning_fragment_count = len(eligible_spanning_fragments)
        observed_fraction = observed_supporting_fragment_fraction(
            supporting_fragment_count,
            spanning_fragment_count,
        )
        filter_reasons = _call_filter_reasons(
            supporting_fragment_count=supporting_fragment_count,
            spanning_fragment_count=spanning_fragment_count,
            observed_fraction=observed_fraction,
            partial_observation=representative.is_partial_observation,
            forward_support_count=forward_support_count,
            reverse_support_count=reverse_support_count,
            discordant_fragment_count=len(consensus.discordant_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            supporting_fragment_count=supporting_fragment_count,
            spanning_fragment_count=spanning_fragment_count,
            observed_supporting_fragment_fraction=observed_fraction,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
            forward_support_count=forward_support_count,
            reverse_support_count=reverse_support_count,
            concordant_fragment_count=len(consensus.concordant_fragment_ids),
            single_mate_fragment_count=len(consensus.single_mate_fragment_ids),
            discordant_fragment_count=len(consensus.discordant_fragment_ids),
            unresolved_fragment_count=len(consensus.unresolved_fragment_ids),
        )
        calls.append(call)
        representatives.extend(
            _sorted_representatives(representative_map[key].values())
        )

    calls.sort(key=_sort_key)
    representatives.sort(key=_representative_sort_key)
    return calls, representatives


def _itd_call_key(itd: ITD) -> ITDCallKey:
    return (
        itd.insertion.start,
        itd.tandem_start,
        itd.tandem_sequence,
        itd.spacer_prefix,
        itd.spacer_suffix,
        itd.insertion.trailing,
    )


def _representative_itd(itds: list[ITD]) -> ITD:
    return itds[0]


def _sort_key(call: ITDCall) -> tuple[int, int, str, str, str]:
    return (
        call.itd.insertion.start,
        call.itd.tandem_start,
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
    right_end = min(len(reference), itd.tandem_start + flank_size)
    right = "".join(
        observed_bases.get(position, "-")
        for position in range(itd.tandem_start, right_end)
    )
    return f"{left}[{itd.tandem_sequence}]{right}"


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
    min_tandem_length: int,
    require_in_frame: bool,
    evidence_filter: InsertionEvidenceFilter | None,
) -> tuple[
    dict[ITDCallKey, list[ITD]],
    SupportRepresentativeMap,
    dict[ITDCallKey, FragmentConsensusSupport],
]:
    alignments = list(alignments)
    grouped_itds: dict[ITDCallKey, list[ITD]] = defaultdict(list)
    representative_map: dict[ITDCallKey, dict[str, tuple[ITD, Alignment]]] = (
        defaultdict(dict)
    )
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
            evidence_filter=evidence_filter,
        )
        for insertion in insertions:
            itd = classify_exact_itd(
                insertion,
                reference,
                min_tandem_length=min_tandem_length,
            )
            if itd is None:
                continue
            observations.append(_CandidateObservation(itd=itd, alignment=alignment))

    consensus_by_key = _fragment_consensus_support(observations, alignments)
    for observation in observations:
        itd = observation.itd
        alignment = observation.alignment

        key = _itd_call_key(itd)
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
        )
        if alignment.fragment_id in consensus_by_key[key].supporting_fragment_ids:
            fragment_ids_by_signature[key][signature].add(alignment.fragment_id)
            insert_sequences_by_fragment[key][alignment.fragment_id][
                itd.insertion.sequence
            ] = _sequence_mismatches(itd.insertion.sequence, expected_sequence)

    finalized_map: SupportRepresentativeMap = defaultdict(dict)
    for key, alignments_by_signature in representative_map.items():
        insert_sequence_supports = _insert_sequence_supports(
            insert_sequences_by_fragment[key]
        )
        for signature, (itd, alignment) in alignments_by_signature.items():
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
            )

    return grouped_itds, finalized_map, consensus_by_key


def _collect_fuzzy_itd_support(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_mismatches: int,
    min_insert_length: int,
    min_tandem_length: int,
    require_in_frame: bool,
    evidence_filter: InsertionEvidenceFilter | None,
) -> tuple[
    dict[ITDCallKey, list[ITD]],
    SupportRepresentativeMap,
    dict[ITDCallKey, FragmentConsensusSupport],
]:
    alignments = list(alignments)
    grouped_itds: dict[ITDCallKey, list[ITD]] = defaultdict(list)
    representative_map: dict[ITDCallKey, dict[str, tuple[ITD, Alignment]]] = (
        defaultdict(dict)
    )
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
            evidence_filter=evidence_filter,
        )
        for insertion in insertions:
            itd = classify_fuzzy_itd(
                insertion,
                reference,
                max_mismatches=max_mismatches,
                min_tandem_length=min_tandem_length,
            )
            if itd is None:
                continue

            observations.append(_CandidateObservation(itd=itd, alignment=alignment))

    consensus_by_key = _fragment_consensus_support(observations, alignments)
    for observation in observations:
        itd = observation.itd
        alignment = observation.alignment
        insertion = itd.insertion

        key = _itd_call_key(itd)
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
        )
        is_supporting = (
            alignment.fragment_id in consensus_by_key[key].supporting_fragment_ids
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
        for signature, (itd, alignment) in alignments_by_signature.items():
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
            )

    return grouped_itds, finalized_map, consensus_by_key


def _fragment_consensus_support(
    observations: Iterable[_CandidateObservation],
    alignments: Iterable[Alignment],
) -> dict[ITDCallKey, FragmentConsensusSupport]:
    """Reconcile R1/R2 evidence before assigning fragment support.

    Concordant fragments have the same candidate and observed inserted sequence
    in both directions.
    Single-mate fragments have one candidate and no opposite mate spanning the
    breakpoint. Discordant fragments have incompatible candidates or a
    breakpoint-spanning opposite mate without the candidate. Multiple candidate
    observations in one direction, or different inserted sequences assigned to
    one candidate across mates, are unresolved. Discordant and unresolved
    fragments are ineligible for both support and local coverage.
    """
    observations = list(observations)
    alignments = list(alignments)
    candidates_by_fragment_site_direction: dict[
        str, dict[int, dict[str, set[tuple[ITDCallKey, str]]]]
    ] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    all_keys: set[ITDCallKey] = set()
    for observation in observations:
        key = _itd_call_key(observation.itd)
        all_keys.add(key)
        candidates_by_fragment_site_direction[observation.alignment.fragment_id][
            observation.itd.insertion.start
        ][observation.alignment.direction].add(
            (key, observation.itd.insertion.sequence)
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
            "concordant": set(),
            "single_mate": set(),
            "discordant": set(),
            "unresolved": set(),
        }
        for key in all_keys
    }

    for fragment_id, sites in candidates_by_fragment_site_direction.items():
        for site, candidates_by_direction in sites.items():
            forward_candidates = candidates_by_direction.get("forward", set())
            reverse_candidates = candidates_by_direction.get("reverse", set())
            candidate_keys = {
                key for key, _ in forward_candidates | reverse_candidates
            }
            if len(forward_candidates) > 1 or len(reverse_candidates) > 1:
                for key in candidate_keys:
                    category_sets[key]["unresolved"].add(fragment_id)
                continue

            if forward_candidates and reverse_candidates:
                forward_key, forward_sequence = next(iter(forward_candidates))
                reverse_key, reverse_sequence = next(iter(reverse_candidates))
                if forward_key != reverse_key:
                    for key in candidate_keys:
                        category_sets[key]["discordant"].add(fragment_id)
                elif forward_sequence == reverse_sequence:
                    key = forward_key
                    category_sets[key]["supporting"].add(fragment_id)
                    category_sets[key]["concordant"].add(fragment_id)
                else:
                    category_sets[forward_key]["unresolved"].add(fragment_id)
                continue

            key = next(iter(candidate_keys))
            candidate_direction = "forward" if forward_candidates else "reverse"
            opposite_direction = (
                "reverse" if candidate_direction == "forward" else "forward"
            )
            opposite_spans_site = any(
                spans_insertion_site(alignment, site)
                for alignment in alignments_by_fragment_direction[fragment_id].get(
                    opposite_direction,
                    [],
                )
            )
            if opposite_spans_site:
                category_sets[key]["discordant"].add(fragment_id)
            else:
                category_sets[key]["supporting"].add(fragment_id)
                category_sets[key]["single_mate"].add(fragment_id)

    return {
        key: FragmentConsensusSupport(
            supporting_fragment_ids=frozenset(categories["supporting"]),
            concordant_fragment_ids=frozenset(categories["concordant"]),
            single_mate_fragment_ids=frozenset(categories["single_mate"]),
            discordant_fragment_ids=frozenset(categories["discordant"]),
            unresolved_fragment_ids=frozenset(categories["unresolved"]),
        )
        for key, categories in category_sets.items()
    }


def _set_best_representative(
    representatives: dict[str, tuple[ITD, Alignment]],
    signature: str,
    itd: ITD,
    alignment: Alignment,
) -> None:
    current = representatives.get(signature)
    candidate_key = (
        _sequence_mismatches(
            itd.insertion.sequence,
            _expected_insertion_sequence(itd),
        ),
        alignment.read_id,
    )
    if current is None:
        representatives[signature] = (itd, alignment)
        return

    current_itd, current_alignment = current
    current_key = (
        _sequence_mismatches(
            current_itd.insertion.sequence,
            _expected_insertion_sequence(current_itd),
        ),
        current_alignment.read_id,
    )
    if candidate_key < current_key:
        representatives[signature] = (itd, alignment)


def _insert_sequence_supports(
    sequences_by_fragment: dict[str, dict[str, int]],
) -> tuple[InsertSequenceSupport, ...]:
    fragment_ids_by_sequence: dict[str, set[str]] = defaultdict(set)
    mismatches_by_sequence: dict[str, int] = {}

    for fragment_id, sequence_mismatches in sequences_by_fragment.items():
        if len(sequence_mismatches) != 1:
            raise ValueError(
                "supporting fragment has unresolved inserted sequences: "
                f"{fragment_id}"
            )
        sequence, mismatches = next(iter(sequence_mismatches.items()))
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
        representative.itd.tandem_start,
        representative.itd.tandem_sequence,
        representative.itd.spacer_prefix,
        representative.itd.spacer_suffix,
        -representative.support_count,
        representative.signature,
        representative.alignment.read_id,
    )


def _call_filter_reasons(
    *,
    supporting_fragment_count: int,
    spanning_fragment_count: int,
    observed_fraction: float,
    partial_observation: bool,
    forward_support_count: int,
    reverse_support_count: int,
    discordant_fragment_count: int,
    unresolved_fragment_count: int,
    filters: ITDFilter,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if partial_observation:
        reasons.append("PARTIAL_OBSERVATION")
    if supporting_fragment_count == 0 and discordant_fragment_count:
        reasons.append("ONLY_DISCORDANT_MATE_EVIDENCE")
    if supporting_fragment_count == 0 and unresolved_fragment_count:
        reasons.append("ONLY_UNRESOLVED_MATE_EVIDENCE")
    if supporting_fragment_count < filters.min_supporting_fragment_count:
        reasons.append("LOW_SUPPORT")
    if spanning_fragment_count < filters.min_spanning_fragment_count:
        reasons.append("LOW_COVERAGE")
    if observed_fraction < filters.min_observed_supporting_fragment_fraction:
        reasons.append("LOW_SUPPORTING_FRAGMENT_FRACTION")
    directional_observations = forward_support_count + reverse_support_count
    if directional_observations >= filters.min_directional_observations:
        largest_direction_fraction = max(
            forward_support_count,
            reverse_support_count,
        ) / directional_observations
        if largest_direction_fraction > filters.max_single_direction_fraction:
            reasons.append("DIRECTION_BIAS")
    return tuple(reasons)


def _direction_support_counts(
    itds: Iterable[ITD],
    supporting_fragment_ids: frozenset[str],
) -> tuple[int, int]:
    fragments_by_direction: dict[str, set[str]] = {
        "forward": set(),
        "reverse": set(),
    }
    for itd in itds:
        if itd.insertion.fragment_id not in supporting_fragment_ids:
            continue
        fragments_by_direction[itd.insertion.direction].add(
            itd.insertion.fragment_id
        )
    return (
        len(fragments_by_direction["forward"]),
        len(fragments_by_direction["reverse"]),
    )


def _expected_insertion_sequence(itd: ITD) -> str:
    return f"{itd.spacer_prefix}{itd.tandem_sequence}{itd.spacer_suffix}"


def _resolved_min_tandem_length(
    min_insert_length: int,
    min_tandem_length: int | None,
) -> int:
    if min_insert_length < 1:
        raise ValueError("min_insert_length must be at least 1")
    resolved = min_insert_length if min_tandem_length is None else min_tandem_length
    if resolved < 1:
        raise ValueError("min_tandem_length must be at least 1")
    return resolved


def _warn_deprecated_call_attribute(old_name: str, new_name: str) -> None:
    scientific_warning = (
        " The fraction is not a validated VAF or allelic ratio."
        if old_name in {"vaf", "min_vaf"}
        else ""
    )
    warnings.warn(
        f"{old_name} is deprecated; use {new_name}.{scientific_warning}",
        DeprecationWarning,
        stacklevel=3,
    )
