import pytest

from itdiscover.alleles import CanonicalInsertionAllele
from itdiscover.calls import ITDCall, ITDFilter
from itdiscover.insertions import Insertion
from itdiscover.itds import ITD
from itdiscover.pipeline import (
    call_exact_itds_from_fragments,
    call_fuzzy_itds_from_fragments,
)
from itdiscover.reads import Fragment, SequencingRead


def make_read(
    read_id: str,
    fragment_id: str,
    sequence: str,
    direction: str,
    quality: int = 40,
) -> SequencingRead:
    return SequencingRead(
        read_id=read_id,
        fragment_id=fragment_id,
        sequence=sequence,
        qualities=(quality,) * len(sequence),
        direction=direction,
    )


def make_fragment(
    fragment_id: str,
    forward_sequence: str,
    reverse_sequence: str,
    forward_quality: int = 40,
    reverse_quality: int = 40,
) -> Fragment:
    return Fragment(
        fragment_id=fragment_id,
        forward_read=make_read(
            f"{fragment_id}/1",
            fragment_id,
            forward_sequence,
            "forward",
            quality=forward_quality,
        ),
        reverse_read=make_read(
            f"{fragment_id}/2",
            fragment_id,
            reverse_sequence,
            "reverse",
            quality=reverse_quality,
        ),
    )


def test_call_exact_itds_from_fragments_reports_observed_fragment_fraction() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            f"itd-fragment-{index}",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
        )
        for index in range(1, 4)
    ] + [
        make_fragment(f"wt-fragment-{index}", reference, reference)
        for index in range(1, 8)
    ]

    assert call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    ) == [
        ITDCall(
            itd=ITD(
                insertion=Insertion(
                    read_id="itd-fragment-1/1",
                    fragment_id="itd-fragment-1",
                    start=2,
                    sequence="CCCGGG",
                    direction="forward",
                ),
                copied_segment_start=3,
                copied_segment_sequence="CCCGGG",
            ),
            canonical_allele=CanonicalInsertionAllele(
                start=2,
                sequence="CCCGGG",
            ),
            mutant_fragment_count=3,
            wild_type_fragment_count=7,
        )
    ]


def test_call_exact_itds_from_fragments_filters_low_quality_reads() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            "low-quality-itd",
            "AAACCCGGGCCCGGGTTT",
            reference,
            forward_quality=10,
        ),
        make_fragment("wt-fragment", reference, reference),
    ]

    assert (
        call_exact_itds_from_fragments(
            fragments,
            reference,
            min_read_length=12,
            min_mean_quality=30,
        )
        == []
    )


def test_call_exact_itds_from_fragments_keeps_passing_mate_when_other_mate_fails() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            "reverse-itd",
            reference,
            "AAACCCGGGCCCGGGTTT",
            forward_quality=10,
        ),
        make_fragment("wt-fragment", reference, reference),
    ]

    calls = call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    )

    assert len(calls) == 1
    assert calls[0].itd.insertion.read_id == "reverse-itd/2"
    assert calls[0].itd.insertion.direction == "reverse"
    assert calls[0].mutant_fragment_count == 1
    assert calls[0].informative_fragment_count == 2
    assert calls[0].observed_mutant_fragment_fraction == 0.5


def test_call_exact_itds_from_fragments_excludes_failed_mate_from_support_but_keeps_passing_mate_for_coverage() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            "failed-itd-passing-wt",
            "AAACCCGGGCCCGGGTTT",
            reference,
            forward_quality=10,
        ),
        make_fragment(
            "passing-itd",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
        ),
    ]

    calls = call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    )

    assert len(calls) == 1
    assert calls[0].mutant_fragment_count == 1
    assert calls[0].informative_fragment_count == 2
    assert calls[0].observed_mutant_fragment_fraction == 0.5


def test_call_exact_itds_from_fragments_trims_terminal_ns_before_calling() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            "itd-fragment",
            "NAAACCCGGGCCCGGGTTTN",
            "AAACCCGGGCCCGGGTTT",
        ),
        make_fragment("wt-fragment", reference, reference),
    ]

    calls = call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    )

    assert len(calls) == 1
    assert calls[0].mutant_fragment_count == 1
    assert calls[0].informative_fragment_count == 2
    assert calls[0].observed_mutant_fragment_fraction == 0.5


def test_call_exact_itds_from_fragments_counts_overlapping_mates_once() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            "itd-fragment",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
        ),
        make_fragment("wt-fragment", reference, reference),
    ]

    calls = call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    )

    assert len(calls) == 1
    assert calls[0].mutant_fragment_count == 1
    assert calls[0].informative_fragment_count == 2
    assert calls[0].observed_mutant_fragment_fraction == 0.5


def test_call_exact_itds_from_fragments_rejects_lowercase_reference() -> None:
    with pytest.raises(ValueError, match="reference contains invalid bases"):
        call_exact_itds_from_fragments([], "AAAccc")


def test_call_fuzzy_itds_from_fragments_reports_observed_fragment_fraction() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            f"itd-fragment-{index}",
            "AAACCCGGGCCCGGATTT",
            "AAACCCGGGCCCGGATTT",
        )
        for index in range(1, 4)
    ] + [
        make_fragment(f"wt-fragment-{index}", reference, reference)
        for index in range(1, 8)
    ]

    assert call_fuzzy_itds_from_fragments(
        fragments,
        reference,
        max_copy_mismatch_rate=1 / 6,
        min_read_length=12,
        min_mean_quality=30,
    ) == [
        ITDCall(
            itd=ITD(
                insertion=Insertion(
                    read_id="itd-fragment-1/1",
                    fragment_id="itd-fragment-1",
                    start=8,
                    sequence="CCCGGA",
                    direction="forward",
                ),
                copied_segment_start=3,
                copied_segment_sequence="CCCGGG",
            ),
            canonical_allele=CanonicalInsertionAllele(
                start=2,
                sequence="CCCGGG",
            ),
            mutant_fragment_count=3,
            wild_type_fragment_count=7,
        )
    ]


def test_call_fuzzy_itds_from_fragments_rejects_itds_over_threshold() -> None:
    reference = "AAACCCGGGTTT"
    fragments = [
        make_fragment(
            "itd-fragment",
            "AAACCCGGGCCCGGATTT",
            reference,
        ),
        make_fragment("wt-fragment", reference, reference),
    ]

    assert (
        call_fuzzy_itds_from_fragments(
            fragments,
            reference,
            max_copy_mismatch_rate=0,
            min_read_length=12,
            min_mean_quality=30,
        )
        == []
    )


def test_fragment_pipeline_propagates_short_out_of_frame_event_policy() -> None:
    reference = "GGGATGCCCTACTTT"
    mutant = "GGGATGCCCACCCTACTTT"
    fragments = [make_fragment("mutant", mutant, reference)]

    assert call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
        min_insert_length=4,
        min_copied_segment_length=3,
    ) == []
    calls = call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=12,
        min_mean_quality=30,
        min_insert_length=4,
        min_copied_segment_length=3,
        require_in_frame=False,
    )

    assert len(calls) == 1
    assert calls[0].itd.copied_segment_sequence == "CCC"


def test_call_fuzzy_itds_from_fragments_rejects_invalid_mismatch_rate() -> None:
    with pytest.raises(
        ValueError,
        match="max_copy_mismatch_rate must be between 0 and 1",
    ):
        call_fuzzy_itds_from_fragments(
            [],
            "AAACCCGGGTTT",
            max_copy_mismatch_rate=-0.01,
        )
