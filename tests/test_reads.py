import pytest

from itdiscover.reads import (
    Fragment,
    ReadTrimSettings,
    SequencingRead,
    orient_read,
    passes_read_filters,
    preprocess_fragments,
    preprocess_fragments_with_metrics,
    preprocess_reads,
    trim_primers,
    trim_terminal_ns,
)
from itdiscover.sequences import reverse_complement

_SequencingRead = SequencingRead


def SequencingRead(*args, **kwargs):
    if args:
        read_id, sequence, qualities, direction = args
        kwargs = {
            "read_id": read_id,
            "fragment_id": kwargs.pop("fragment_id", read_id),
            "sequence": sequence,
            "qualities": qualities,
            "direction": direction,
            **kwargs,
        }
    else:
        kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _SequencingRead(**kwargs)


def test_sequencing_read_requires_quality_for_each_base() -> None:
    with pytest.raises(ValueError, match="equal length"):
        SequencingRead(
            read_id="read-1",
            sequence="ACGT",
            qualities=(30, 30, 30),
            direction="forward",
        )


def test_sequencing_read_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        SequencingRead(
            read_id="read-1",
            sequence="ACgT",
            qualities=(30, 30, 30, 30),
            direction="forward",
        )


def test_reverse_complement_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        reverse_complement("ACgT")


def test_orient_read_keeps_forward_reads_as_is() -> None:
    read = orient_read(
        read_id="r1",
        sequence="ACGTN",
        qualities=(10, 20, 30, 40, 35),
        direction="forward",
        fragment_id="r1",
    )

    assert read == SequencingRead(
        read_id="r1",
        sequence="ACGTN",
        qualities=(10, 20, 30, 40, 35),
        direction="forward",
    )


def test_orient_read_rejects_lowercase_forward_reads() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        orient_read(
            read_id="r1",
            sequence="ACgTN",
            qualities=(10, 20, 30, 40, 35),
            direction="forward",
            fragment_id="r1",
        )


def test_orient_read_reverse_complements_reverse_reads() -> None:
    read = orient_read(
        read_id="r2",
        sequence="ACGTT",
        qualities=(10, 20, 30, 40, 35),
        direction="reverse",
        fragment_id="r2",
    )

    assert read == SequencingRead(
        read_id="r2",
        sequence="AACGT",
        qualities=(35, 40, 30, 20, 10),
        direction="reverse",
    )


def test_orient_read_rejects_lowercase_reverse_reads() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        orient_read(
            read_id="r2",
            sequence="ACgTT",
            qualities=(10, 20, 30, 40, 35),
            direction="reverse",
            fragment_id="r2",
        )


def test_trim_terminal_ns_preserves_internal_ns_and_quality_alignment() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="NNACNGTN",
        qualities=(5, 6, 30, 31, 32, 33, 34, 7),
        direction="forward",
    )

    assert trim_terminal_ns(read) == SequencingRead(
        read_id="read-1",
        sequence="ACNGT",
        qualities=(30, 31, 32, 33, 34),
        direction="forward",
    )


def test_trim_primers_trims_forward_reads() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="GGGCCCGGGTTT",
        qualities=tuple(range(12)),
        direction="forward",
    )

    assert trim_primers(
        read,
        ReadTrimSettings(
            forward_primer="CCC",
        ),
    ) == SequencingRead(
        read_id="read-1",
        sequence="GGGTTT",
        qualities=(6, 7, 8, 9, 10, 11),
        direction="forward",
    )


def test_trim_primers_trims_raw_r2_reverse_primer_after_orientation() -> None:
    raw_r2 = "AAACCCAAACCCGGGTTT"
    read = orient_read(
        read_id="read-2",
        fragment_id="read-2",
        sequence=raw_r2,
        qualities=tuple(range(len(raw_r2))),
        direction="reverse",
    )

    assert trim_primers(
        read,
        ReadTrimSettings(reverse_primer="AAACCC"),
    ) == SequencingRead(
        read_id="read-2",
        sequence="AAACCCGGGTTT",
        qualities=tuple(range(17, 5, -1)),
        direction="reverse",
    )


def test_trim_primers_discards_reads_without_primer() -> None:
    read = SequencingRead(
        read_id="read-3",
        sequence="GGGTTT",
        qualities=(10, 11, 12, 13, 14, 15),
        direction="forward",
    )

    assert trim_primers(
        read,
        ReadTrimSettings(forward_primer="CCC"),
    ) is None


def test_passes_read_filters_uses_trimmed_length_and_mean_quality() -> None:
    passing = SequencingRead(
        read_id="passing",
        sequence="ACGT",
        qualities=(30, 31, 32, 33),
        direction="forward",
    )
    low_quality = SequencingRead(
        read_id="low-quality",
        sequence="ACGT",
        qualities=(10, 10, 10, 10),
        direction="forward",
    )

    assert passes_read_filters(passing, min_length=4, min_mean_quality=30)
    assert not passes_read_filters(passing, min_length=5, min_mean_quality=30)
    assert not passes_read_filters(low_quality, min_length=4, min_mean_quality=30)


def test_fragment_requires_matching_fragment_ids_and_directions() -> None:
    with pytest.raises(ValueError, match="reverse read fragment_id"):
        Fragment(
            fragment_id="fragment-1",
            forward_read=SequencingRead("fragment-1/1", "ACGT", (30,) * 4, "forward", fragment_id="fragment-1"),
            reverse_read=SequencingRead("fragment-2/2", "ACGT", (30,) * 4, "reverse", fragment_id="fragment-2"),
        )


def test_preprocess_reads_keeps_distinct_fragments_separate() -> None:
    reads = [
        SequencingRead("read-1", "ACGT", (30, 30, 30, 30), "forward"),
        SequencingRead("read-2", "ACGT", (31, 31, 31, 31), "forward"),
        SequencingRead("read-3", "ACGT", (32, 32, 32, 32), "reverse"),
    ]

    assert preprocess_reads(reads, min_length=4, min_mean_quality=30) == reads


def test_preprocess_reads_trims_filters_and_collapses() -> None:
    reads = [
        SequencingRead("read-1", "NACGTN", (5, 30, 30, 30, 30, 5), "forward"),
        SequencingRead("read-2", "ACGT", (31, 31, 31, 31), "forward"),
        SequencingRead("read-3", "NNNN", (40, 40, 40, 40), "forward"),
        SequencingRead("read-4", "ACGT", (10, 10, 10, 10), "forward"),
    ]

    assert preprocess_reads(reads, min_length=4, min_mean_quality=30) == [
        SequencingRead("read-1", "ACGT", (30, 30, 30, 30), "forward"),
        SequencingRead("read-2", "ACGT", (31, 31, 31, 31), "forward"),
    ]


def test_preprocess_fragments_returns_passing_reads_from_paired_fragments() -> None:
    fragment = Fragment(
        fragment_id="fragment-1",
        forward_read=SequencingRead("fragment-1/1", "ACGT", (30,) * 4, "forward", fragment_id="fragment-1"),
        reverse_read=SequencingRead("fragment-1/2", "NNNN", (30,) * 4, "reverse", fragment_id="fragment-1"),
    )

    assert preprocess_fragments([fragment], min_length=4, min_mean_quality=30) == [
        SequencingRead("fragment-1/1", "ACGT", (30,) * 4, "forward", fragment_id="fragment-1")
    ]


def test_preprocess_fragments_retains_rejection_and_direction_metrics() -> None:
    fragments = [
        Fragment(
            fragment_id="passing",
            forward_read=SequencingRead(
                "passing/1",
                "AAAACGT",
                (40,) * 7,
                "forward",
                fragment_id="passing",
            ),
            reverse_read=SequencingRead(
                "passing/2",
                "ACGT",
                (40,) * 4,
                "reverse",
                fragment_id="passing",
            ),
        ),
        Fragment(
            fragment_id="missing-primer",
            forward_read=SequencingRead(
                "missing-primer/1",
                "TTTACGT",
                (40,) * 7,
                "forward",
                fragment_id="missing-primer",
            ),
            reverse_read=SequencingRead(
                "missing-primer/2",
                "ACGT",
                (40,) * 4,
                "reverse",
                fragment_id="missing-primer",
            ),
        ),
    ]

    result = preprocess_fragments_with_metrics(
        fragments,
        min_length=4,
        min_mean_quality=30,
        trimming=ReadTrimSettings(forward_primer="AAA"),
    )

    assert result.metrics.input_fragment_count == 2
    assert result.metrics.input_read_count == 4
    assert result.metrics.primer_retained_forward_reads == 1
    assert result.metrics.primer_retained_reverse_reads is None
    assert result.metrics.primer_failed_read_count == 1
    assert result.metrics.passing_read_count == 3
    assert result.metrics.passing_forward_read_count == 1
    assert result.metrics.passing_reverse_read_count == 2
    assert result.metrics.usable_fragment_count == 2
