import pytest

from itdiscover.insertions import (
    Alignment,
    Insertion,
    InsertionEvidenceFilter,
    extract_insertions,
)

_Alignment = Alignment
_Insertion = Insertion


def Alignment(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Alignment(**kwargs)


def Insertion(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Insertion(**kwargs)


def test_alignment_requires_equal_aligned_lengths() -> None:
    with pytest.raises(ValueError, match="equal length"):
        Alignment(
            read_id="read-1",
            read_sequence="ACGT",
            aligned_read="ACGT",
            aligned_reference="ACG",
            direction="forward",
        )


def test_alignment_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="read_sequence contains invalid bases"):
        Alignment(
            read_id="read-1",
            read_sequence="ACgT",
            aligned_read="ACGT",
            aligned_reference="ACGT",
            direction="forward",
        )


def test_alignment_rejects_invalid_alignment_characters() -> None:
    with pytest.raises(ValueError, match="aligned_read contains invalid bases"):
        Alignment(
            read_id="read-1",
            read_sequence="ACGT",
            aligned_read="ACXT",
            aligned_reference="ACGT",
            direction="forward",
        )


def test_insertion_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="sequence contains invalid bases"):
        Insertion(
            read_id="read-1",
            start=2,
            sequence="ACgT",
            direction="forward",
        )


def test_extracts_internal_in_frame_insertion() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCCCCCGGGTTT",
        aligned_read="AAACCCCCCGGGTTT",
        aligned_reference="AAA------GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCCCC",
            direction="forward",
            trailing=False,
        )
    ]


def test_extracts_multiple_insertions_from_one_alignment() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCCCCCGGGTTTTTTAAA",
        aligned_read="AAACCCCCCGGGTTTTTTAAA",
        aligned_reference="AAA------GGG------AAA",
        direction="forward",
    )

    assert extract_insertions(alignment) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCCCC",
            direction="forward",
        ),
        Insertion(
            read_id="read-1",
            start=5,
            sequence="TTTTTT",
            direction="forward",
        ),
    ]


def test_filters_short_insertions_by_default() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCGGGTTT",
        aligned_read="AAACCGGGTTT",
        aligned_reference="AAA--GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment) == []
    assert extract_insertions(alignment, min_length=2, require_in_frame=False) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CC",
            direction="forward",
        )
    ]


def test_filters_internal_out_of_frame_insertions() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCCCGGGTTT",
        aligned_read="AAACCCCGGGTTT",
        aligned_reference="AAA----GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment, min_length=4) == []
    assert extract_insertions(alignment, min_length=4, require_in_frame=False) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCC",
            direction="forward",
        )
    ]


def test_allows_out_of_frame_trailing_insertions() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="CCCCAAA",
        aligned_read="CCCCAAA",
        aligned_reference="----AAA",
        direction="reverse",
    )

    assert extract_insertions(alignment, min_length=4) == [
        Insertion(
            read_id="read-1",
            start=-1,
            sequence="CCCC",
            direction="reverse",
            trailing=True,
        )
    ]


def test_filters_insertions_with_ambiguous_bases() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACNCGGGTTT",
        aligned_read="AAACNCGGGTTT",
        aligned_reference="AAA---GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment, min_length=3) == []


def test_allows_isolated_q29_insert_base_when_insert_mean_passes() -> None:
    alignment = Alignment(
        read_id="isolated-q29",
        read_sequence="AAACCCCCCGGG",
        aligned_read="AAACCCCCCGGG",
        aligned_reference="AAA------GGG",
        direction="forward",
        aligned_qualities=(40, 40, 40, 29, 40, 40, 40, 40, 40, 40, 40, 40),
    )

    assert len(
        extract_insertions(
            alignment,
            evidence_filter=InsertionEvidenceFilter(),
        )
    ) == 1


def test_filters_insertions_with_low_quality_junction_anchor() -> None:
    alignment = Alignment(
        read_id="low-anchor-quality",
        read_sequence="AAACCCCCCGGG",
        aligned_read="AAACCCCCCGGG",
        aligned_reference="AAA------GGG",
        direction="forward",
        aligned_qualities=(40, 40, 29, 40, 40, 40, 40, 40, 40, 40, 40, 40),
    )

    assert extract_insertions(
        alignment,
        evidence_filter=InsertionEvidenceFilter(),
    ) == []


def test_filters_insertions_with_low_mean_insert_quality() -> None:
    alignment = Alignment(
        read_id="low-insert-mean",
        read_sequence="AAACCCCCCGGG",
        aligned_read="AAACCCCCCGGG",
        aligned_reference="AAA------GGG",
        direction="forward",
        aligned_qualities=(40, 40, 40, 29, 29, 29, 29, 29, 29, 40, 40, 40),
    )

    assert extract_insertions(
        alignment,
        evidence_filter=InsertionEvidenceFilter(),
    ) == []


def test_filters_insertions_with_one_very_low_quality_insert_base() -> None:
    alignment = Alignment(
        read_id="low-insert-base",
        read_sequence="AAACCCCCCGGG",
        aligned_read="AAACCCCCCGGG",
        aligned_reference="AAA------GGG",
        direction="forward",
        aligned_qualities=(40, 40, 40, 14, 40, 40, 40, 40, 40, 40, 40, 40),
    )

    assert extract_insertions(
        alignment,
        evidence_filter=InsertionEvidenceFilter(),
    ) == []
