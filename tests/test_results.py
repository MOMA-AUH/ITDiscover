from dataclasses import replace

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
            tandem_start=3,
            tandem_sequence="CCCGGG",
            orientation="downstream",
        ),
        supporting_fragment_count=1,
        spanning_fragment_count=10,
        observed_supporting_fragment_fraction=0.1,
        status="FAIL",
        filter_reasons=("LOW_SUPPORT",),
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


def test_analysis_error_is_not_reported_as_negative() -> None:
    result = error_sample_result("broken-sample", ValueError("invalid FASTQ"))

    assert result.analysis_status == "error"
    assert result.qc_status == "fail"
    assert result.outcome == "indeterminate"
    assert result.qc_reasons == ("ANALYSIS_ERROR",)
    assert result.error_message == "invalid FASTQ"
