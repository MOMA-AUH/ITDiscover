import pytest

from itdiscover.calls import (
    ITDCall,
    ITDFilter,
    call_exact_itds,
    call_fuzzy_itds_with_representatives,
)
from itdiscover.insertions import Alignment, Insertion
from itdiscover.itds import ITD


def make_alignment(
    read_id: str,
    read_sequence: str,
    aligned_read: str,
    aligned_reference: str,
    direction: str = "forward",
    fragment_id: str | None = None,
) -> Alignment:
    return Alignment(
        read_id=read_id,
        fragment_id=fragment_id or read_id,
        read_sequence=read_sequence,
        aligned_read=aligned_read,
        aligned_reference=aligned_reference,
        direction=direction,
    )


def permissive_filters() -> ITDFilter:
    return ITDFilter(
        min_supporting_fragment_count=1,
        min_spanning_fragment_count=0,
        min_observed_supporting_fragment_fraction=0,
    )


def make_insertion(
    read_id: str,
    start: int,
    sequence: str,
    direction: str = "forward",
) -> Insertion:
    return Insertion(
        read_id=read_id,
        fragment_id=read_id,
        start=start,
        sequence=sequence,
        direction=direction,
    )


def test_itd_call_rejects_fraction_inconsistent_with_fragment_counts() -> None:
    itd = ITD(
        insertion=make_insertion("read", start=2, sequence="CCCGGG"),
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
    )

    with pytest.raises(ValueError, match="must equal"):
        ITDCall(
            itd=itd,
            supporting_fragment_count=1,
            spanning_fragment_count=2,
            observed_supporting_fragment_fraction=0.25,
        )
    with pytest.raises(ValueError, match="must equal"):
        ITDCall(
            itd=itd,
            supporting_fragment_count=1,
            spanning_fragment_count=2,
            observed_supporting_fragment_fraction=float("nan"),
        )


def test_call_exact_itds_reports_supporting_and_spanning_fragment_fraction() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            f"itd-read-{index}",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
        )
        for index in range(1, 4)
    ] + [
        make_alignment(
            f"wt-read-{index}",
            reference,
            reference,
            reference,
            direction="reverse",
        )
        for index in range(1, 8)
    ]

    assert call_exact_itds(alignments, reference) == [
        ITDCall(
            itd=ITD(
                insertion=make_insertion(
                    "itd-read-1",
                    start=8,
                    sequence="CCCGGG",
                ),
                tandem_start=3,
                tandem_sequence="CCCGGG",
                orientation="upstream",
            ),
            supporting_fragment_count=3,
            spanning_fragment_count=10,
            observed_supporting_fragment_fraction=0.3,
        )
    ]


def test_call_exact_itds_counts_concordant_and_single_mate_fragments() -> None:
    reference = "AAACCCGGGTTT"
    mutant = "AAACCCGGGCCCGGGTTT"
    aligned_reference = "AAACCCGGG------TTT"
    alignments = [
        make_alignment(
            "concordant/1",
            mutant,
            mutant,
            aligned_reference,
            fragment_id="concordant",
        ),
        make_alignment(
            "concordant/2",
            mutant,
            mutant,
            aligned_reference,
            direction="reverse",
            fragment_id="concordant",
        ),
        make_alignment(
            "single/1",
            mutant,
            mutant,
            aligned_reference,
            fragment_id="single",
        ),
    ]

    call = call_exact_itds(alignments, reference, filters=permissive_filters())[0]

    assert call.supporting_fragment_count == 2
    assert call.spanning_fragment_count == 2
    assert call.observed_supporting_fragment_fraction == 1
    assert call.concordant_fragment_count == 1
    assert call.single_mate_fragment_count == 1
    assert call.discordant_fragment_count == 0
    assert call.unresolved_fragment_count == 0


def test_call_exact_itds_excludes_mutant_wt_mate_disagreement() -> None:
    reference = "AAACCCGGGTTT"
    mutant = "AAACCCGGGCCCGGGTTT"
    alignments = [
        make_alignment(
            "discordant/1",
            mutant,
            mutant,
            "AAACCCGGG------TTT",
            fragment_id="discordant",
        ),
        make_alignment(
            "discordant/2",
            reference,
            reference,
            reference,
            direction="reverse",
            fragment_id="discordant",
        ),
    ]

    call = call_exact_itds(alignments, reference, filters=permissive_filters())[0]

    assert call.supporting_fragment_count == 0
    assert call.spanning_fragment_count == 0
    assert call.observed_supporting_fragment_fraction == 0
    assert call.discordant_fragment_count == 1
    assert call.filter_reasons == (
        "ONLY_DISCORDANT_MATE_EVIDENCE",
        "LOW_SUPPORT",
    )


def test_call_exact_itds_excludes_mutually_incompatible_mate_candidates() -> None:
    reference = "AAACCCGGGTTTAAACCCGGG"
    aligned_reference = "AAACCCGGG------TTTAAACCCGGG"
    alignments = [
        make_alignment(
            "discordant/1",
            "AAACCCGGGCCCGGGTTTAAACCCGGG",
            "AAACCCGGGCCCGGGTTTAAACCCGGG",
            aligned_reference,
            fragment_id="discordant",
        ),
        make_alignment(
            "discordant/2",
            "AAACCCGGGTTTAAATTTAAACCCGGG",
            "AAACCCGGGTTTAAATTTAAACCCGGG",
            aligned_reference,
            direction="reverse",
            fragment_id="discordant",
        ),
    ]

    calls = call_exact_itds(alignments, reference, filters=permissive_filters())

    assert len(calls) == 2
    assert {call.itd.tandem_sequence for call in calls} == {"CCCGGG", "TTTAAA"}
    assert all(call.supporting_fragment_count == 0 for call in calls)
    assert all(call.spanning_fragment_count == 0 for call in calls)
    assert all(call.discordant_fragment_count == 1 for call in calls)
    assert all(
        call.filter_reasons
        == ("ONLY_DISCORDANT_MATE_EVIDENCE", "LOW_SUPPORT")
        for call in calls
    )


def test_call_exact_itds_marks_multiple_same_mate_candidates_unresolved() -> None:
    reference = "AAACCCGGGTTTAAACCCGGG"
    aligned_reference = "AAACCCGGG------TTTAAACCCGGG"
    alignments = [
        make_alignment(
            "ambiguous/alignment-1",
            "AAACCCGGGCCCGGGTTTAAACCCGGG",
            "AAACCCGGGCCCGGGTTTAAACCCGGG",
            aligned_reference,
            fragment_id="ambiguous",
        ),
        make_alignment(
            "ambiguous/alignment-2",
            "AAACCCGGGTTTAAATTTAAACCCGGG",
            "AAACCCGGGTTTAAATTTAAACCCGGG",
            aligned_reference,
            fragment_id="ambiguous",
        ),
    ]

    calls = call_exact_itds(alignments, reference, filters=permissive_filters())

    assert len(calls) == 2
    assert all(call.supporting_fragment_count == 0 for call in calls)
    assert all(call.spanning_fragment_count == 0 for call in calls)
    assert all(call.unresolved_fragment_count == 1 for call in calls)
    assert all(
        call.filter_reasons
        == ("ONLY_UNRESOLVED_MATE_EVIDENCE", "LOW_SUPPORT")
        for call in calls
    )


def test_call_fuzzy_itds_does_not_hide_mate_sequence_disagreement() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            "fuzzy-disagreement/1",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
            fragment_id="fuzzy-disagreement",
        ),
        make_alignment(
            "fuzzy-disagreement/2",
            "AAACCCGGGCCCGGATTT",
            "AAACCCGGGCCCGGATTT",
            "AAACCCGGG------TTT",
            direction="reverse",
            fragment_id="fuzzy-disagreement",
        ),
    ]

    calls, representatives = call_fuzzy_itds_with_representatives(
        alignments,
        reference,
        max_mismatches=1,
        filters=permissive_filters(),
    )

    assert len(calls) == 1
    assert calls[0].supporting_fragment_count == 0
    assert calls[0].spanning_fragment_count == 0
    assert calls[0].unresolved_fragment_count == 1
    assert calls[0].filter_reasons == (
        "ONLY_UNRESOLVED_MATE_EVIDENCE",
        "LOW_SUPPORT",
    )
    assert all(
        representative.support_count == 0 for representative in representatives
    )
    assert all(
        representative.insert_sequence_supports == ()
        for representative in representatives
    )


def test_call_exact_itds_keeps_distinct_breakpoints_separate() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            f"upstream-representation-{index}",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
        )
        for index in range(1, 3)
    ] + [
        make_alignment(
            f"downstream-representation-{index}",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAA------CCCGGGTTT",
            direction="reverse",
        )
        for index in range(1, 4)
    ] + [
        make_alignment(
            f"wt-read-{index}",
            reference,
            reference,
            reference,
        )
        for index in range(1, 8)
    ]

    calls = call_exact_itds(alignments, reference)

    assert len(calls) == 2
    assert [call.itd.insertion.start for call in calls] == [2, 8]
    assert [call.supporting_fragment_count for call in calls] == [3, 2]
    assert [call.spanning_fragment_count for call in calls] == [12, 12]
    assert [call.observed_supporting_fragment_fraction for call in calls] == [
        3 / 12,
        2 / 12,
    ]


def test_call_exact_itds_are_independent_of_input_order_for_distinct_breakpoints() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            "upstream-breakpoint",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
        ),
        make_alignment(
            "downstream-breakpoint",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAA------CCCGGGTTT",
        ),
    ]

    calls = call_exact_itds(alignments, reference)
    reversed_calls = call_exact_itds(list(reversed(alignments)), reference)

    assert calls == reversed_calls
    assert [
        (
            call.itd.insertion.start,
            call.supporting_fragment_count,
            call.observed_supporting_fragment_fraction,
        )
        for call in calls
    ] == [(2, 1, 0.5), (8, 1, 0.5)]


def test_call_exact_itds_counts_overlapping_mates_once_per_fragment() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="fragment-1/1",
            fragment_id="fragment-1",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAA------CCCGGGTTT",
            direction="forward",
        ),
        Alignment(
            read_id="fragment-1/2",
            fragment_id="fragment-1",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAA------CCCGGGTTT",
            direction="reverse",
        ),
        Alignment(
            read_id="fragment-2/1",
            fragment_id="fragment-2",
            read_sequence=reference,
            aligned_read=reference,
            aligned_reference=reference,
            direction="forward",
        ),
    ]

    calls = call_exact_itds(alignments, reference)

    assert len(calls) == 1
    assert calls[0].supporting_fragment_count == 1
    assert calls[0].spanning_fragment_count == 2
    assert calls[0].observed_supporting_fragment_fraction == 0.5


def test_call_exact_itds_reports_unique_supporting_sequences() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="fragment-1/1",
            fragment_id="fragment-1",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAACCCGGG------TTT",
            direction="forward",
        ),
        Alignment(
            read_id="fragment-2/1",
            fragment_id="fragment-2",
            read_sequence="GGGAAACCCGGGCCCGGGTTT",
            aligned_read="GGGAAACCCGGGCCCGGGTTT",
            aligned_reference="---AAACCCGGG------TTT",
            direction="forward",
        ),
        Alignment(
            read_id="fragment-3/1",
            fragment_id="fragment-3",
            read_sequence="AAACTTGGGCCCGGGTTT",
            aligned_read="AAACTTGGGCCCGGGTTT",
            aligned_reference="AAACCCGGG------TTT",
            direction="forward",
        ),
    ] + [
        make_alignment(f"wt-read-{index}", reference, reference, reference)
        for index in range(1, 4)
    ]

    calls = call_exact_itds(alignments, reference)

    assert len(calls) == 1
    assert calls[0].supporting_fragment_count == 3
    assert calls[0].spanning_fragment_count == 6
    assert calls[0].observed_supporting_fragment_fraction == 0.5


def test_call_exact_itds_ignores_non_itd_insertions() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            f"insertion-read-{index}",
            "AAACCCGGGAAAAAATTT",
            "AAACCCGGGAAAAAATTT",
            "AAACCCGGG------TTT",
        )
        for index in range(1, 5)
    ] + [
        make_alignment(f"wt-read-{index}", reference, reference, reference)
        for index in range(1, 7)
    ]

    assert call_exact_itds(alignments, reference) == []


def test_call_exact_itds_ignores_reference_match_away_from_breakpoint() -> None:
    reference = "AAACCCGGGTTTAAATTTGGG"
    alignments = [
        make_alignment(
            "remote-copy-read",
            "AAACCCGGGTTTAAACCCGGGTTTGGG",
            "AAACCCGGGTTTAAACCCGGGTTTGGG",
            "AAACCCGGGTTTAAA------TTTGGG",
        )
    ]

    assert call_exact_itds(alignments, reference) == []


def test_call_exact_itds_marks_read_edge_itds_as_partial_observations() -> None:
    reference = "AAACCCGGG"
    alignment = make_alignment(
        "partial-itd",
        "AAACCCGGGCCCGGG",
        "AAACCCGGGCCCGGG",
        "AAACCCGGG------",
    )

    calls = call_exact_itds(
        [alignment],
        reference,
        filters=ITDFilter(
            min_supporting_fragment_count=1,
            min_spanning_fragment_count=0,
            min_observed_supporting_fragment_fraction=0,
        ),
    )

    assert len(calls) == 1
    assert calls[0].itd.tandem_sequence == "CCCGGG"
    assert calls[0].itd.is_partial_observation
    assert calls[0].status == "FAIL"
    assert calls[0].filter_reasons == ("PARTIAL_OBSERVATION",)


def test_call_exact_itds_reconstructs_fully_observed_long_itd() -> None:
    reference = "AAACCCGGGTTTAAACCCGGGTTT"
    tandem = "CCCGGGTTTAAA"
    alignment = make_alignment(
        "long-itd",
        f"AAA{tandem}{reference[3:]}",
        f"AAA{tandem}{reference[3:]}",
        f"AAA{'-' * len(tandem)}{reference[3:]}",
    )

    calls = call_exact_itds(
        [alignment],
        reference,
        filters=ITDFilter(
            min_supporting_fragment_count=1,
            min_spanning_fragment_count=0,
            min_observed_supporting_fragment_fraction=0,
        ),
    )

    assert len(calls) == 1
    assert calls[0].itd.tandem_sequence == tandem
    assert calls[0].itd.length == len(tandem)
    assert not calls[0].itd.is_partial_observation
    assert calls[0].status == "PASS"


def test_call_exact_itds_returns_sorted_calls() -> None:
    reference = "AAACCCGGGTTTAAACCC"
    alignments = [
            make_alignment(
                "later-itd",
                "AAACCCGGGTTTAAACCCAAACCC",
                "AAACCCGGGTTTAAACCCAAACCC",
                "AAACCCGGGTTTAAACCC------",
            ),
        make_alignment(
            "earlier-itd",
            "AAACCCGGGCCCGGGTTTAAACCC",
            "AAACCCGGGCCCGGGTTTAAACCC",
            "AAACCCGGG------TTTAAACCC",
        ),
    ]

    calls = call_exact_itds(alignments, reference, min_insert_length=3)

    assert [call.itd.insertion.read_id for call in calls] == [
        "earlier-itd",
        "later-itd",
    ]


def test_lower_insert_threshold_also_allows_three_base_tandem_by_default() -> None:
    reference = "AAACCCGGGTTT"
    alignment = make_alignment(
        "three-base-itd",
        "AAACCCCCCGGGTTT",
        "AAACCCCCCGGGTTT",
        "AAA---CCCGGGTTT",
    )

    calls = call_exact_itds(
        [alignment],
        reference,
        min_insert_length=3,
    )

    assert len(calls) == 1
    assert calls[0].itd.tandem_sequence == "CCC"
    assert call_exact_itds(
        [alignment],
        reference,
        min_insert_length=3,
        min_tandem_length=6,
    ) == []


def test_out_of_frame_insertion_filter_is_explicitly_configurable() -> None:
    reference = "AAACCCGGGTTT"
    alignment = make_alignment(
        "out-of-frame-itd",
        "AAACCCACCCGGGTTT",
        "AAACCCACCCGGGTTT",
        "AAA----CCCGGGTTT",
    )

    assert call_exact_itds(
        [alignment],
        reference,
        min_insert_length=4,
        min_tandem_length=3,
    ) == []
    calls = call_exact_itds(
        [alignment],
        reference,
        min_insert_length=4,
        min_tandem_length=3,
        require_in_frame=False,
    )

    assert len(calls) == 1
    assert calls[0].itd.tandem_sequence == "CCC"
    assert calls[0].itd.spacer_suffix == "A"


def test_call_exact_itds_marks_call_as_fail_when_support_threshold_is_not_met() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            "itd-read",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
        ),
        make_alignment("wt-read", reference, reference, reference),
    ]

    calls = call_exact_itds(
        alignments,
        reference,
        filters=ITDFilter(
            min_supporting_fragment_count=2,
            min_spanning_fragment_count=0,
            min_observed_supporting_fragment_fraction=0,
        ),
    )

    assert len(calls) == 1
    assert calls[0].status == "FAIL"
    assert calls[0].filter_reasons == ("LOW_SUPPORT",)


def test_default_call_filters_do_not_pass_a_single_fragment_candidate() -> None:
    reference = "AAACCCGGGTTT"
    alignment = make_alignment(
        "single-fragment",
        "AAACCCGGGCCCGGGTTT",
        "AAACCCGGGCCCGGGTTT",
        "AAACCCGGG------TTT",
    )

    calls = call_exact_itds([alignment], reference)

    assert calls[0].status == "FAIL"
    assert calls[0].filter_reasons == ("LOW_SUPPORT", "LOW_COVERAGE")


def test_filters_on_observed_supporting_fragment_fraction() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            "supporting-fragment",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
        )
    ] + [
        make_alignment(f"spanning-fragment-{index}", reference, reference, reference)
        for index in range(9)
    ]

    calls = call_exact_itds(
        alignments,
        reference,
        filters=ITDFilter(
            min_supporting_fragment_count=1,
            min_spanning_fragment_count=10,
            min_observed_supporting_fragment_fraction=0.2,
        ),
    )

    assert calls[0].observed_supporting_fragment_fraction == 0.1
    assert calls[0].filter_reasons == ("LOW_SUPPORTING_FRAGMENT_FRACTION",)


def test_call_exact_itds_marks_one_direction_only_support_as_direction_biased() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        make_alignment(
            f"forward-only-{index}",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGGCCCGGGTTT",
            "AAACCCGGG------TTT",
        )
        for index in range(5)
    ]

    calls = call_exact_itds(
        alignments,
        reference,
        filters=ITDFilter(
            min_supporting_fragment_count=1,
            min_spanning_fragment_count=0,
            min_observed_supporting_fragment_fraction=0,
        ),
    )

    assert calls[0].forward_support_count == 5
    assert calls[0].reverse_support_count == 0
    assert calls[0].status == "FAIL"
    assert calls[0].filter_reasons == ("DIRECTION_BIAS",)


def test_call_fuzzy_itds_reports_exact_and_fuzzy_only_support_counts() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="exact-fragment/1",
            fragment_id="exact-fragment",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAA------CCCGGGTTT",
            direction="forward",
        ),
        Alignment(
            read_id="fuzzy-fragment/1",
            fragment_id="fuzzy-fragment",
            read_sequence="AAACCCGGACCCGGGTTT",
            aligned_read="AAACCCGGACCCGGGTTT",
            aligned_reference="AAA------CCCGGGTTT",
            direction="forward",
        ),
    ]

    calls, representatives = call_fuzzy_itds_with_representatives(
        alignments,
        reference,
        max_mismatches=1,
    )

    assert len(calls) == 1
    assert calls[0].supporting_fragment_count == 2
    assert len(representatives) == 1
    assert representatives[0].support_count == 2
    assert representatives[0].exact_support_count == 1
    assert representatives[0].fuzzy_only_support_count == 1
    assert representatives[0].fuzzy_example_sequence == "CCCGGA"
    assert representatives[0].mismatches == 0
    assert representatives[0].insert_sequence_supports[0].sequence == "CCCGGG"
    assert representatives[0].insert_sequence_supports[0].support_count == 1
    assert representatives[0].insert_sequence_supports[0].mismatches == 0
    assert representatives[0].insert_sequence_supports[1].sequence == "CCCGGA"
    assert representatives[0].insert_sequence_supports[1].support_count == 1
    assert representatives[0].insert_sequence_supports[1].mismatches == 1


def test_call_fuzzy_itds_groups_spacer_itd_reads_with_copied_segment_mismatches() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="exact-fragment-1/1",
            fragment_id="exact-fragment-1",
            read_sequence="AAATTACCCGGGACTCCCGGGTTT",
            aligned_read="AAATTACCCGGGACTCCCGGGTTT",
            aligned_reference="AAA------------CCCGGGTTT",
            direction="forward",
        ),
        Alignment(
            read_id="exact-fragment-2/1",
            fragment_id="exact-fragment-2",
            read_sequence="AAATTACCCGGGACTCCCGGGTTT",
            aligned_read="AAATTACCCGGGACTCCCGGGTTT",
            aligned_reference="AAA------------CCCGGGTTT",
            direction="forward",
        ),
        Alignment(
            read_id="fuzzy-fragment/1",
            fragment_id="fuzzy-fragment",
            read_sequence="AAATTACCCGGAACTCCCGGGTTT",
            aligned_read="AAATTACCCGGAACTCCCGGGTTT",
            aligned_reference="AAA------------CCCGGGTTT",
            direction="forward",
        ),
    ]

    calls, representatives = call_fuzzy_itds_with_representatives(
        alignments,
        reference,
        max_mismatches=1,
    )

    assert len(calls) == 1
    assert calls[0].supporting_fragment_count == 3
    assert calls[0].itd == ITD(
        insertion=Insertion(
            read_id="exact-fragment-1/1",
            fragment_id="exact-fragment-1",
            start=2,
            sequence="TTACCCGGGACT",
            direction="forward",
        ),
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
        spacer_prefix="TTA",
        spacer_suffix="ACT",
    )
    assert len(representatives) == 1
    assert representatives[0].support_count == 3
    assert representatives[0].exact_support_count == 2
    assert representatives[0].fuzzy_only_support_count == 1
    assert representatives[0].fuzzy_example_sequence == "TTACCCGGAACT"
    assert representatives[0].insert_sequence_supports[0].sequence == "TTACCCGGGACT"
    assert representatives[0].insert_sequence_supports[0].support_count == 2
    assert representatives[0].insert_sequence_supports[0].mismatches == 0
    assert representatives[0].insert_sequence_supports[1].sequence == "TTACCCGGAACT"
    assert representatives[0].insert_sequence_supports[1].support_count == 1
    assert representatives[0].insert_sequence_supports[1].mismatches == 1
