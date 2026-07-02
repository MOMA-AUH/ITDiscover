"""Build reportable ITD calls from aligned reads."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from .coverage import interbase_coverage, variant_allele_frequency
from .insertions import Alignment, Insertion, extract_insertions
from .itds import ITD, classify_exact_itd, classify_fuzzy_itd


@dataclass(frozen=True)
class ITDCall:
    """An ITD with support, coverage, and allele frequency."""

    itd: ITD
    support_count: int
    coverage: int
    vaf: float
    status: str = "PASS"
    filter_reasons: tuple[str, ...] = ()

    @property
    def passes_filters(self) -> bool:
        """Return whether the call passes the configured thresholds."""
        return self.status == "PASS"


@dataclass(frozen=True)
class ITDFilter:
    """Thresholds used to label exact-match ITD calls."""

    min_support_count: int = 1
    min_coverage: int = 0
    min_vaf: float = 0.0

    def __post_init__(self) -> None:
        if self.min_support_count < 1:
            raise ValueError("min_support_count must be at least 1")
        if self.min_coverage < 0:
            raise ValueError("min_coverage must not be negative")
        if self.min_vaf < 0:
            raise ValueError("min_vaf must not be negative")


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


ITDCallKey = tuple[int, str, str, str, bool]
SupportRepresentativeMap = dict[ITDCallKey, dict[str, UniqueSupportRepresentative]]


def call_exact_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
    filters: ITDFilter = ITDFilter(),
) -> list[ITDCall]:
    """Call exact-match ITDs and attach support, coverage, and VAF."""
    calls, _ = call_exact_itds_with_representatives(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        filters=filters,
    )
    return calls


def call_fuzzy_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_mismatches: int,
    min_insert_length: int = 6,
    filters: ITDFilter = ITDFilter(),
) -> list[ITDCall]:
    """Call fuzzy-match ITDs and attach support, coverage, and VAF."""
    if max_mismatches < 0:
        raise ValueError("max_mismatches must not be negative")
    calls, _ = call_fuzzy_itds_with_representatives(
        alignments,
        reference,
        max_mismatches=max_mismatches,
        min_insert_length=min_insert_length,
        filters=filters,
    )
    return calls


def call_fuzzy_itds_with_representatives(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_mismatches: int,
    min_insert_length: int = 6,
    filters: ITDFilter = ITDFilter(),
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call fuzzy-match ITDs and retain one alignment per unique support pattern."""
    if max_mismatches < 0:
        raise ValueError("max_mismatches must not be negative")
    alignments = list(alignments)
    coverage_by_site = interbase_coverage(alignments)
    grouped_itds, representative_map = _collect_fuzzy_itd_support(
        alignments,
        reference,
        max_mismatches=max_mismatches,
        min_insert_length=min_insert_length,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for itds in grouped_itds.values():
        representative = _representative_itd(itds)
        key = _itd_call_key(representative)
        support_count = len({itd.insertion.fragment_id for itd in itds})
        coverage = coverage_by_site.get(representative.insertion.start, 0)
        vaf = variant_allele_frequency(support_count, coverage)
        filter_reasons = _call_filter_reasons(
            support_count=support_count,
            coverage=coverage,
            vaf=vaf,
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            support_count=support_count,
            coverage=coverage,
            vaf=vaf,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
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
    filters: ITDFilter = ITDFilter(),
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call exact-match ITDs and retain one alignment per unique support pattern."""
    alignments = list(alignments)
    coverage_by_site = interbase_coverage(alignments)
    grouped_itds, representative_map = _collect_exact_itd_support(
        alignments,
        reference,
        min_insert_length=min_insert_length,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for itds in grouped_itds.values():
        representative = _representative_itd(itds)
        key = _itd_call_key(representative)
        support_count = len({itd.insertion.fragment_id for itd in itds})
        coverage = coverage_by_site.get(representative.insertion.start, 0)
        vaf = variant_allele_frequency(support_count, coverage)
        filter_reasons = _call_filter_reasons(
            support_count=support_count,
            coverage=coverage,
            vaf=vaf,
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            support_count=support_count,
            coverage=coverage,
            vaf=vaf,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
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
) -> tuple[dict[ITDCallKey, list[ITD]], SupportRepresentativeMap]:
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

    for alignment in alignments:
        insertions = extract_insertions(alignment, min_length=min_insert_length)
        for insertion in insertions:
            itd = classify_exact_itd(insertion, reference)
            if itd is None:
                continue

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

    return grouped_itds, finalized_map


def _collect_fuzzy_itd_support(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    max_mismatches: int,
    min_insert_length: int,
) -> tuple[dict[ITDCallKey, list[ITD]], SupportRepresentativeMap]:
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

    for alignment in alignments:
        insertions = extract_insertions(alignment, min_length=min_insert_length)
        for insertion in insertions:
            itd = classify_fuzzy_itd(
                insertion,
                reference,
                max_mismatches=max_mismatches,
            )
            if itd is None:
                continue

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
            fragment_ids_by_signature[key][signature].add(alignment.fragment_id)
            insert_sequences_by_fragment[key][alignment.fragment_id][
                insertion.sequence
            ] = sequence_mismatches

            if sequence_mismatches == 0:
                exact_fragment_ids_by_signature[key][signature].add(
                    alignment.fragment_id
                )
                continue

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

    return grouped_itds, finalized_map


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
    support_count: int,
    coverage: int,
    vaf: float,
    filters: ITDFilter,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if support_count < filters.min_support_count:
        reasons.append("LOW_SUPPORT")
    if coverage < filters.min_coverage:
        reasons.append("LOW_COVERAGE")
    if vaf < filters.min_vaf:
        reasons.append("LOW_VAF")
    return tuple(reasons)


def _expected_insertion_sequence(itd: ITD) -> str:
    return f"{itd.spacer_prefix}{itd.tandem_sequence}{itd.spacer_suffix}"
