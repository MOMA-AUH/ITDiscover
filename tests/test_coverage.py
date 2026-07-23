import pytest

from itdiscover.coverage import (
    covered_reference_positions,
    interbase_coverage,
    interbase_fragment_ids,
    observed_mutant_fragment_fraction,
    spans_insertion_site,
)
from itdiscover.insertions import Alignment

_Alignment = Alignment


def Alignment(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Alignment(**kwargs)


def test_covered_reference_positions_excludes_insertions_and_deletions() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AACCCT",
        aligned_read="AA---C-T",
        aligned_reference="AAGGG-CT",
        direction="forward",
    )

    assert covered_reference_positions(alignment) == {0, 1, 6}


def test_internal_site_requires_coverage_on_both_sides() -> None:
    spanning = Alignment(
        read_id="spanning",
        read_sequence="AAACCC",
        aligned_read="AAACCC",
        aligned_reference="AAACCC",
        direction="forward",
    )
    deleted_right_side = Alignment(
        read_id="deleted",
        read_sequence="AAACC",
        aligned_read="AAA-CC",
        aligned_reference="AAACCC",
        direction="forward",
    )

    assert spans_insertion_site(spanning, 2)
    assert not spans_insertion_site(deleted_right_side, 2)


def test_edge_sites_require_terminal_reference_base_coverage() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AACCC",
        aligned_read="-AACCC",
        aligned_reference="AAACCC",
        direction="forward",
    )

    assert not spans_insertion_site(alignment, -1)
    assert spans_insertion_site(alignment, 5)


def test_interbase_coverage_counts_distinct_fragments() -> None:
    alignments = [
        Alignment(
            read_id=f"full-read-{index}",
            read_sequence="AAACCC",
            aligned_read="AAACCC",
            aligned_reference="AAACCC",
            direction="forward",
        )
        for index in range(1, 4)
    ] + [
        Alignment(
            read_id=f"partial-read-{index}",
            read_sequence="CCC",
            aligned_read="---CCC",
            aligned_reference="AAACCC",
            direction="reverse",
        )
        for index in range(1, 3)
    ]

    assert interbase_coverage(alignments) == {
        -1: 3,
        0: 3,
        1: 3,
        2: 3,
        3: 5,
        4: 5,
        5: 5,
    }


def test_interbase_coverage_counts_overlapping_mates_once_per_fragment() -> None:
    alignments = [
        Alignment(
            read_id="fragment-1/1",
            fragment_id="fragment-1",
            read_sequence="AAACCC",
            aligned_read="AAACCC",
            aligned_reference="AAACCC",
            direction="forward",
        ),
        Alignment(
            read_id="fragment-1/2",
            fragment_id="fragment-1",
            read_sequence="AAACCC",
            aligned_read="AAACCC",
            aligned_reference="AAACCC",
            direction="reverse",
        ),
        Alignment(
            read_id="fragment-2/1",
            fragment_id="fragment-2",
            read_sequence="AAACCC",
            aligned_read="AAACCC",
            aligned_reference="AAACCC",
            direction="forward",
        ),
    ]

    assert interbase_coverage(alignments) == {
        -1: 2,
        0: 2,
        1: 2,
        2: 2,
        3: 2,
        4: 2,
        5: 2,
    }
    assert interbase_fragment_ids(alignments)[2] == frozenset(
        {"fragment-1", "fragment-2"}
    )


def test_spans_insertion_site_rejects_out_of_range_site() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAA",
        aligned_read="AAA",
        aligned_reference="AAA",
        direction="forward",
    )

    with pytest.raises(ValueError, match="outside"):
        spans_insertion_site(alignment, 3)


def test_observed_mutant_fragment_fraction_returns_fraction() -> None:
    assert observed_mutant_fragment_fraction(2, 8) == 0.25
    assert observed_mutant_fragment_fraction(0, 8) == 0.0
    assert observed_mutant_fragment_fraction(0, 0) == 0.0


def test_observed_mutant_fragment_fraction_rejects_impossible_counts() -> None:
    with pytest.raises(ValueError, match="must not be negative"):
        observed_mutant_fragment_fraction(-1, 10)
    with pytest.raises(ValueError, match="must not be negative"):
        observed_mutant_fragment_fraction(1, -10)
    with pytest.raises(ValueError, match="must not exceed"):
        observed_mutant_fragment_fraction(11, 10)
    with pytest.raises(ValueError, match="must not exceed"):
        observed_mutant_fragment_fraction(1, 0)
