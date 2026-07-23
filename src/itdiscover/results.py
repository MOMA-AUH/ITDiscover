"""Sample-level analysis, quality-control, and outcome contracts."""

from dataclasses import dataclass
from statistics import median
from typing import Literal

from .calls import ITDCall
from .coverage import interbase_coverage
from .insertions import Alignment
from .reads import PreprocessingMetrics

AnalysisStatus = Literal["complete", "error"]
QCStatus = Literal["pass", "warn", "fail"]
Outcome = Literal[
    "ITD detected",
    "no passing ITD detected",
    "indeterminate",
]


@dataclass(frozen=True)
class SampleQCThresholds:
    """Minimum assay-adequacy thresholds required for a negative outcome."""

    min_usable_fragment_count: int = 10
    min_passing_reads_per_direction: int = 1
    min_alignment_pass_fraction: float = 0.80
    min_median_interbase_coverage: int = 10
    min_primer_retention_fraction: float = 0.80

    def __post_init__(self) -> None:
        if self.min_usable_fragment_count < 1:
            raise ValueError("min_usable_fragment_count must be at least 1")
        if self.min_passing_reads_per_direction < 0:
            raise ValueError("min_passing_reads_per_direction must not be negative")
        if not 0 <= self.min_alignment_pass_fraction <= 1:
            raise ValueError("min_alignment_pass_fraction must be between 0 and 1")
        if self.min_median_interbase_coverage < 0:
            raise ValueError("min_median_interbase_coverage must not be negative")
        if not 0 <= self.min_primer_retention_fraction <= 1:
            raise ValueError("min_primer_retention_fraction must be between 0 and 1")


@dataclass(frozen=True)
class AlignmentMetrics:
    """Counts retained through read-to-reference alignment filtering."""

    attempted_read_count: int
    passing_read_count: int
    passing_forward_read_count: int
    passing_reverse_read_count: int
    passing_fragment_count: int

    @property
    def pass_fraction(self) -> float:
        """Return the fraction of attempted alignments that passed."""
        if self.attempted_read_count == 0:
            return 0.0
        return self.passing_read_count / self.attempted_read_count


@dataclass(frozen=True)
class CoverageMetrics:
    """Distribution of distinct-fragment inter-base coverage across the target."""

    minimum: int
    median: float
    maximum: int


@dataclass(frozen=True)
class SampleResult:
    """Top-level sample result separating execution, QC, and biological outcome."""

    sample_id: str
    analysis_status: AnalysisStatus
    qc_status: QCStatus
    outcome: Outcome
    qc_reasons: tuple[str, ...]
    preprocessing: PreprocessingMetrics | None
    alignment: AlignmentMetrics | None
    coverage: CoverageMetrics | None
    passing_call_count: int
    filtered_candidate_count: int
    error_message: str | None = None


def alignment_metrics(
    attempted_read_count: int,
    passing_alignments: list[Alignment],
) -> AlignmentMetrics:
    """Build alignment-retention metrics from filtered alignments."""
    return AlignmentMetrics(
        attempted_read_count=attempted_read_count,
        passing_read_count=len(passing_alignments),
        passing_forward_read_count=sum(
            alignment.direction == "forward" for alignment in passing_alignments
        ),
        passing_reverse_read_count=sum(
            alignment.direction == "reverse" for alignment in passing_alignments
        ),
        passing_fragment_count=len(
            {alignment.fragment_id for alignment in passing_alignments}
        ),
    )


def coverage_metrics(
    alignments: list[Alignment],
    reference_length: int,
) -> CoverageMetrics:
    """Return min/median/max coverage including uncovered inter-base sites."""
    coverage_by_site = interbase_coverage(alignments)
    depths = [
        coverage_by_site.get(site, 0)
        for site in range(-1, reference_length)
    ]
    if not depths:
        return CoverageMetrics(minimum=0, median=0.0, maximum=0)
    return CoverageMetrics(
        minimum=min(depths),
        median=float(median(depths)),
        maximum=max(depths),
    )


def build_sample_result(
    *,
    sample_id: str,
    calls: list[ITDCall],
    preprocessing: PreprocessingMetrics,
    alignment: AlignmentMetrics,
    coverage: CoverageMetrics,
    thresholds: SampleQCThresholds,
) -> SampleResult:
    """Evaluate QC and derive an outcome from passing calls and assay adequacy."""
    reasons: list[str] = []
    if preprocessing.input_fragment_count == 0:
        reasons.append("NO_INPUT_FRAGMENTS")
    if preprocessing.usable_fragment_count < thresholds.min_usable_fragment_count:
        reasons.append("LOW_USABLE_FRAGMENT_COUNT")
    if (
        alignment.passing_forward_read_count
        < thresholds.min_passing_reads_per_direction
    ):
        reasons.append("LOW_FORWARD_READ_COUNT")
    if (
        alignment.passing_reverse_read_count
        < thresholds.min_passing_reads_per_direction
    ):
        reasons.append("LOW_REVERSE_READ_COUNT")
    if alignment.pass_fraction < thresholds.min_alignment_pass_fraction:
        reasons.append("LOW_ALIGNMENT_PASS_FRACTION")
    if coverage.median < thresholds.min_median_interbase_coverage:
        reasons.append("LOW_MEDIAN_INTERBASE_COVERAGE")

    input_fragments = preprocessing.input_fragment_count
    if input_fragments > 0:
        if preprocessing.primer_retained_forward_reads is not None and (
            preprocessing.primer_retained_forward_reads / input_fragments
            < thresholds.min_primer_retention_fraction
        ):
            reasons.append("LOW_FORWARD_PRIMER_RETENTION")
        if preprocessing.primer_retained_reverse_reads is not None and (
            preprocessing.primer_retained_reverse_reads / input_fragments
            < thresholds.min_primer_retention_fraction
        ):
            reasons.append("LOW_REVERSE_PRIMER_RETENTION")

    passing_call_count = sum(call.passes_filters for call in calls)
    filtered_candidate_count = len(calls) - passing_call_count
    failed_sample_qc = bool(reasons)
    ambiguous_evidence_dominates = any(
        "AMBIGUOUS_EVIDENCE_DOMINATES" in call.filter_reasons
        for call in calls
    )
    if ambiguous_evidence_dominates:
        reasons.append("AMBIGUOUS_ITD_EVIDENCE_DOMINATES")

    if failed_sample_qc or (
        ambiguous_evidence_dominates and passing_call_count == 0
    ):
        qc_status: QCStatus = "fail"
    elif filtered_candidate_count:
        qc_status = "warn"
        reasons.append("FILTERED_ITD_CANDIDATES_PRESENT")
    else:
        qc_status = "pass"
    if qc_status == "fail":
        outcome: Outcome = "indeterminate"
    elif passing_call_count:
        outcome = "ITD detected"
    else:
        outcome = "no passing ITD detected"

    return SampleResult(
        sample_id=sample_id,
        analysis_status="complete",
        qc_status=qc_status,
        outcome=outcome,
        qc_reasons=tuple(reasons),
        preprocessing=preprocessing,
        alignment=alignment,
        coverage=coverage,
        passing_call_count=passing_call_count,
        filtered_candidate_count=filtered_candidate_count,
    )


def error_sample_result(sample_id: str, error: Exception) -> SampleResult:
    """Return an explicit analysis-error result without assay metrics."""
    return SampleResult(
        sample_id=sample_id,
        analysis_status="error",
        qc_status="fail",
        outcome="indeterminate",
        qc_reasons=("ANALYSIS_ERROR",),
        preprocessing=None,
        alignment=None,
        coverage=None,
        passing_call_count=0,
        filtered_candidate_count=0,
        error_message=str(error),
    )
