"""Inter-base coverage and observed supporting-fragment fractions."""

import warnings

from collections import defaultdict
from collections.abc import Iterable

from .insertions import Alignment


def covered_reference_positions(alignment: Alignment) -> set[int]:
    """Return zero-based reference positions covered by read bases."""
    covered: set[int] = set()
    ref_pos = -1

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        if ref_base == "-":
            continue
        ref_pos += 1
        if read_base != "-":
            covered.add(ref_pos)

    return covered


def reference_length(alignment: Alignment) -> int:
    """Return the ungapped reference length represented by an alignment."""
    return sum(1 for base in alignment.aligned_reference if base != "-")


def spans_insertion_site(alignment: Alignment, site: int) -> bool:
    """Return whether an alignment spans an inter-base insertion site.

    `site` uses the same convention as `Insertion.start`: an insertion at
    `site` occurs after reference base `site`. Leading insertions use `-1`.
    Internal sites require read bases on both sides of the inter-base position.
    """
    ref_length = reference_length(alignment)
    if site < -1 or site >= ref_length:
        raise ValueError("site is outside the aligned reference")

    covered = covered_reference_positions(alignment)
    if site == -1:
        return 0 in covered
    if site == ref_length - 1:
        return site in covered
    return site in covered and site + 1 in covered


def interbase_coverage(alignments: Iterable[Alignment]) -> dict[int, int]:
    """Return fragment-level coverage for every spanned insertion site."""
    return {
        site: len(fragment_ids)
        for site, fragment_ids in interbase_fragment_ids(alignments).items()
    }


def interbase_fragment_ids(
    alignments: Iterable[Alignment],
) -> dict[int, frozenset[str]]:
    """Return distinct fragment IDs spanning each inter-base insertion site."""
    coverage: defaultdict[int, set[str]] = defaultdict(set)

    for alignment in alignments:
        ref_length = reference_length(alignment)
        for site in range(-1, ref_length):
            if spans_insertion_site(alignment, site):
                coverage[site].add(alignment.fragment_id)

    return {
        site: frozenset(fragment_ids)
        for site, fragment_ids in coverage.items()
    }


def observed_supporting_fragment_fraction(
    supporting_fragment_count: int,
    spanning_fragment_count: int,
) -> float:
    """Return supporting fragments divided by locally spanning fragments.

    Counts are post-filter, fragment-deduplicated counts. Zero spanning
    coverage is defined as 0.0 only when there is also no
    supporting evidence.  Positive support without spanning coverage is an
    impossible count combination and must not be represented as a fraction.
    """
    if supporting_fragment_count < 0:
        raise ValueError("supporting_fragment_count must not be negative")
    if spanning_fragment_count < 0:
        raise ValueError("spanning_fragment_count must not be negative")
    if supporting_fragment_count > spanning_fragment_count:
        raise ValueError(
            "supporting_fragment_count must not exceed spanning_fragment_count"
        )
    if spanning_fragment_count == 0:
        return 0.0
    return supporting_fragment_count / spanning_fragment_count


def variant_allele_frequency(supporting_count: int, spanning_count: int) -> float:
    """Deprecated compatibility alias for the observed fragment fraction."""
    warnings.warn(
        "variant_allele_frequency() is deprecated; use "
        "observed_supporting_fragment_fraction(). The result is not a "
        "validated VAF or allelic ratio.",
        DeprecationWarning,
        stacklevel=2,
    )
    return observed_supporting_fragment_fraction(supporting_count, spanning_count)
