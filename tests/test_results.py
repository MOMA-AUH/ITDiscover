from dataclasses import replace

from itdiscover.alleles import CanonicalInsertionAllele
from itdiscover.calls import ITDCall
from itdiscover.insertions import Insertion
from itdiscover.itds import ITD
from itdiscover.reads import PreprocessingMetrics
from itdiscover.results import (
    AlignmentMetrics,
    CoverageMetrics,
    SampleQCThresholds,
    build_sample_result,
    error_sample_result,
)


def preprocessing_metrics(*, usable_fragments: int = 10) -> PreprocessingMetrics:
    return PreprocessingMetrics(
        input_fragment_count=10,
        input_read_count=20,
        primer_retained_forward_reads=None,
        primer_retained_reverse_reads=None,
        primer_failed_read_count=0,
        length_failed_read_count=0,
        quality_failed_read_count=0,
        passing_read_count=20,
        passing_forward_read_count=10,
        passing_reverse_read_count=10,
        usable_fragment_count=usable_fragments,
    )


def alignment_metrics(*, passing_reads: int = 20) -> AlignmentMetrics:
    return AlignmentMetrics(
        attempted_read_count=20,
        passing_read_count=passing_reads,
        passing_forward_read_count=passing_reads // 2,
        passing_reverse_read_count=passing_reads // 2,
        passing_fragment_count=passing_reads // 2,
    )


def test_adequate_sample_without_calls_has_explicit_negative_outcome() -> None:
    result = build_sample_result(
        sample_id="negative-sample",
        calls=[],
        preprocessing=preprocessing_metrics(),
        alignment=alignment_metrics(),
        coverage=CoverageMetrics(minimum=10, median=10, maximum=10),
        thresholds=SampleQCThresholds(),
    )

    assert result.analysis_status == "complete"
    assert result.qc_status == "pass"
    assert result.outcome == "no passing ITD detected"
    assert result.qc_reasons == ()


def test_inadequate_sample_without_calls_is_indeterminate() -> None:
    result = build_sample_result(
        sample_id="failed-sample",
        calls=[],
        preprocessing=preprocessing_metrics(usable_fragments=2),
        alignment=alignment_metrics(passing_reads=4),
        coverage=CoverageMetrics(minimum=0, median=2, maximum=4),
        thresholds=SampleQCThresholds(),
    )

    assert result.qc_status == "fail"
    assert result.outcome == "indeterminate"
    assert result.qc_reasons == (
        "LOW_USABLE_FRAGMENT_COUNT",
        "LOW_ALIGNMENT_PASS_FRACTION",
        "LOW_MEDIAN_INTERBASE_COVERAGE",
    )


def test_configured_primer_retention_is_an_explicit_qc_gate() -> None:
    preprocessing = replace(
        preprocessing_metrics(),
        primer_retained_forward_reads=7,
        primer_retained_reverse_reads=10,
    )

    result = build_sample_result(
        sample_id="primer-failure",
        calls=[],
        preprocessing=preprocessing,
        alignment=alignment_metrics(),
        coverage=CoverageMetrics(minimum=10, median=10, maximum=10),
        thresholds=SampleQCThresholds(min_primer_retention_fraction=0.8),
    )

    assert result.qc_status == "fail"
    assert result.outcome == "indeterminate"
    assert result.qc_reasons == ("LOW_FORWARD_PRIMER_RETENTION",)


def test_filtered_candidate_produces_qc_warning_without_becoming_detected() -> None:
    candidate = ITDCall(
        itd=ITD(
            insertion=Insertion(
                read_id="candidate/1",
                fragment_id="candidate",
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
        mutant_fragment_count=1,
        wild_type_fragment_count=9,
        status="FAIL",
        filter_reasons=("LOW_MUTANT_FRAGMENT_COUNT",),
    )

    result = build_sample_result(
        sample_id="filtered-candidate",
        calls=[candidate],
        preprocessing=preprocessing_metrics(),
        alignment=alignment_metrics(),
        coverage=CoverageMetrics(minimum=10, median=10, maximum=10),
        thresholds=SampleQCThresholds(),
    )

    assert result.qc_status == "warn"
    assert result.outcome == "no passing ITD detected"
    assert result.qc_reasons == ("FILTERED_ITD_CANDIDATES_PRESENT",)
    assert result.filtered_candidate_count == 1


def test_dominant_ambiguous_evidence_makes_no_call_sample_indeterminate() -> None:
    candidate = ITDCall(
        itd=ITD(
            insertion=Insertion(
                read_id="candidate/1",
                fragment_id="candidate",
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
        conflicting_fragment_count=97,
        status="FAIL",
        filter_reasons=("AMBIGUOUS_EVIDENCE_DOMINATES",),
    )

    result = build_sample_result(
        sample_id="ambiguous-candidate",
        calls=[candidate],
        preprocessing=preprocessing_metrics(),
        alignment=alignment_metrics(),
        coverage=CoverageMetrics(minimum=10, median=10, maximum=10),
        thresholds=SampleQCThresholds(),
    )

    assert result.qc_status == "fail"
    assert result.outcome == "indeterminate"
    assert result.qc_reasons == ("AMBIGUOUS_ITD_EVIDENCE_DOMINATES",)


def test_dominant_ambiguous_candidate_warns_when_another_call_passes() -> None:
    ambiguous = ITDCall(
        itd=ITD(
            insertion=Insertion(
                read_id="ambiguous/1",
                fragment_id="ambiguous",
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
        mutant_fragment_count=1,
        wild_type_fragment_count=9,
        unresolved_fragment_count=11,
        status="FAIL",
        filter_reasons=("AMBIGUOUS_EVIDENCE_DOMINATES",),
    )
    passing = replace(
        ambiguous,
        mutant_fragment_count=5,
        wild_type_fragment_count=5,
        unresolved_fragment_count=0,
        status="PASS",
        filter_reasons=(),
    )

    result = build_sample_result(
        sample_id="passing-with-ambiguous-candidate",
        calls=[passing, ambiguous],
        preprocessing=preprocessing_metrics(),
        alignment=alignment_metrics(),
        coverage=CoverageMetrics(minimum=10, median=10, maximum=10),
        thresholds=SampleQCThresholds(),
    )

    assert result.qc_status == "warn"
    assert result.outcome == "ITD detected"
    assert result.qc_reasons == (
        "AMBIGUOUS_ITD_EVIDENCE_DOMINATES",
        "FILTERED_ITD_CANDIDATES_PRESENT",
    )


def test_analysis_error_is_not_reported_as_negative() -> None:
    result = error_sample_result("broken-sample", ValueError("invalid FASTQ"))

    assert result.analysis_status == "error"
    assert result.qc_status == "fail"
    assert result.outcome == "indeterminate"
    assert result.qc_reasons == ("ANALYSIS_ERROR",)
    assert result.error_message == "invalid FASTQ"
