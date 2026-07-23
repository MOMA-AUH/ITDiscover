from pathlib import Path

import pytest

from itdiscover.insertions import Insertion
from itdiscover.itds import (
    ITD,
    TandemSimilarity,
    classify_exact_itd,
    classify_fuzzy_itd,
    score_tandem_similarity,
)

_Insertion = Insertion


def Insertion(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Insertion(**kwargs)


def test_classifies_copy_before_insertion_without_rewriting_breakpoint() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=8,
        sequence="CCCGGG",
        direction="forward",
    )

    itd = classify_exact_itd(insertion, "AAACCCGGGTTT")

    assert itd is not None
    assert itd == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
    )
    assert itd.copied_segment_location == "before"


def test_classifies_copy_after_insertion() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGG",
        direction="reverse",
    )

    itd = classify_exact_itd(insertion, "AAACCCGGGTTT")

    assert itd is not None
    assert itd == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
    )
    assert itd.copied_segment_location == "after"


def test_reports_copied_segment_from_full_length_flt3_reference() -> None:
    reference_path = (
        Path(__file__).parent / "data" / "synthetic_flt3" / "reference.fa"
    )
    reference = "".join(reference_path.read_text().splitlines()[1:])
    insertion = Insertion(
        read_id="flt3-example",
        start=78,
        sequence="AGAGAATATGAATAT",
        direction="forward",
    )

    itd = classify_exact_itd(insertion, reference)

    assert itd is not None
    assert itd.copied_segment_start == 79
    assert itd.copied_segment_end == 93
    assert itd.copied_segment_sequence == "AGAGAATATGAATAT"
    assert itd.copied_segment_location == "after"


def test_uses_most_five_prime_source_interval_when_both_sides_match() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=5,
        sequence="AAAAAA",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAAAAAAAAAAA") == ITD(
        insertion=insertion,
        tandem_start=0,
        tandem_sequence="AAAAAA",
    )


def test_uses_copy_after_insertion_for_leading_repetitive_itd() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=0,
        sequence="AAAAAA",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAAAAAAAAAAA") == ITD(
        insertion=insertion,
        tandem_start=1,
        tandem_sequence="AAAAAA",
    )


def test_classifies_exact_tandem_with_spacers_on_both_sides() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="NNNCCCGGGNN",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        spacer_prefix="NNN",
        spacer_suffix="NN",
    )


def test_does_not_classify_exact_reference_match_away_from_breakpoint() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=14,
        sequence="CCCGGG",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTTAAATTTGGG") is None


def test_exact_classification_matches_zero_mismatch_fuzzy_classification() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="TTACCCGGGACT",
        direction="forward",
    )
    reference = "AAACCCGGGTTT"

    assert classify_exact_itd(insertion, reference) == classify_fuzzy_itd(
        insertion,
        reference,
        max_mismatches=0,
    )


def test_does_not_classify_sequence_absent_from_reference() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="GGGGGG",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") is None


def test_classify_exact_itd_rejects_lowercase_reference() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGG",
        direction="forward",
    )

    with pytest.raises(ValueError, match="reference contains invalid bases"):
        classify_exact_itd(insertion, "AAAcccGGGTTT")


def test_itd_reports_inclusive_tandem_end_and_length() -> None:
    itd = ITD(
        insertion=Insertion(
            read_id="read-1",
            start=8,
            sequence="CCCGGG",
            direction="forward",
        ),
        tandem_start=3,
        tandem_sequence="CCCGGG",
    )

    assert itd.tandem_end == 8
    assert itd.copied_segment_location == "before"
    assert itd.length == 6
    assert itd.spacer_sequence == ""
    assert itd.spacer_length == 0


def test_itd_derives_location_and_validates_legacy_orientation() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=8,
        sequence="CCCGGG",
        direction="forward",
    )

    with pytest.raises(ValueError, match="orientation is inconsistent"):
        ITD(
            insertion=insertion,
            tandem_start=3,
            tandem_sequence="CCCGGG",
            orientation="downstream",
        )
    with pytest.raises(ValueError, match="not adjacent"):
        ITD(
            insertion=insertion,
            tandem_start=0,
            tandem_sequence="CCC",
        )


def test_scores_exact_tandem_similarity() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGG",
        direction="forward",
    )

    similarity = score_tandem_similarity(insertion, "AAACCCGGGTTT")
    
    assert similarity == TandemSimilarity(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        mismatches=0,
    )


def test_scores_best_tandem_similarity_with_one_mismatch() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    similarity = score_tandem_similarity(insertion, "AAACCCGGGTTT")

    assert similarity == TandemSimilarity(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        mismatches=1,
    )
    assert similarity.matches == 5
    assert similarity.identity == 5 / 6


def test_scores_most_five_prime_tandem_on_tie() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=5,
        sequence="AAA",
        direction="forward",
    )

    similarity = score_tandem_similarity(insertion, "AAAAAA")

    assert similarity == TandemSimilarity(
        insertion=insertion,
        tandem_start=0,
        tandem_sequence="AAA",
        mismatches=0,
    )


def test_classifies_exact_tandem_with_longest_copied_segment_and_spacers() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGGTTTAAAA",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGGTTT",
        spacer_prefix="",
        spacer_suffix="AAAA",
    )


def test_classifies_fuzzy_copy_after_insertion_with_one_mismatch() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=1,
    ) == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
    )


def test_classifies_fuzzy_tandem_with_spacers_when_exact_match_exists() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="TTACCCGGGACT",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=1,
    ) == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        spacer_prefix="TTA",
        spacer_suffix="ACT",
    )


def test_does_not_classify_fuzzy_itd_when_mismatches_exceed_threshold() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=0,
    ) is None


def test_does_not_classify_fuzzy_itd_for_sequences_without_a_close_match() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="NNNNNN",
        direction="forward",
    )

    assert classify_fuzzy_itd(insertion, "AAACCCGGGTTT", max_mismatches=3) is None


def test_uses_most_five_prime_window_when_fuzzy_candidates_tie() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=5,
        sequence="AAAAAT",
        direction="forward",
    )

    assert classify_fuzzy_itd(insertion, "AAAAAAA", max_mismatches=1) == ITD(
        insertion=insertion,
        tandem_start=0,
        tandem_sequence="AAAAAA",
    )


def test_classifies_fuzzy_tandem_with_spacers_and_mismatch_in_copied_segment() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="TTACCCGGAACT",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=1,
    ) == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        spacer_prefix="TTA",
        spacer_suffix="ACT",
    )


def test_classify_fuzzy_itd_rejects_negative_mismatch_threshold() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    with pytest.raises(ValueError, match="max_mismatches must not be negative"):
        classify_fuzzy_itd(insertion, "AAACCCGGGTTT", max_mismatches=-1)


def test_minimum_tandem_length_is_configurable_for_three_base_duplication() -> None:
    insertion = Insertion(
        read_id="three-base-itd",
        start=2,
        sequence="CCC",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") is None
    assert classify_exact_itd(
        insertion,
        "AAACCCGGGTTT",
        min_tandem_length=3,
    ) == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCC",
    )


def test_classification_rejects_invalid_minimum_tandem_length() -> None:
    insertion = Insertion(
        read_id="itd",
        start=2,
        sequence="CCC",
        direction="forward",
    )

    with pytest.raises(ValueError, match="min_tandem_length"):
        classify_exact_itd(
            insertion,
            "AAACCCGGGTTT",
            min_tandem_length=0,
        )
