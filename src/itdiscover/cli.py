"""Command-line interface for ITDiscover."""

import argparse
import csv
import hashlib
import html
import sys
from pathlib import Path
from typing import TextIO

from . import __version__
from .alignment import (
    AlignmentEvidenceFilter,
    align_read_to_reference,
    passes_alignment_evidence_filters,
)
from .calls import (
    ITDCall,
    ITDConsolidationSettings,
    ITDFilter,
    UniqueSupportRepresentative,
    call_exact_itds_with_representatives,
    call_fuzzy_itds_with_representatives,
)
from .fastq import read_paired_fastq
from .insertions import (
    Alignment,
    InsertionEvidenceFilter,
)
from .itds import ITD
from .reads import (
    PrimerOrientationError,
    ReadTrimSettings,
    preprocess_fragments_with_metrics,
    validate_primer_orientations,
)
from .results import (
    SampleQCThresholds,
    SampleResult,
    alignment_metrics,
    build_sample_result,
    coverage_metrics,
    error_sample_result,
)


COORDINATE_CONVENTION = (
    "Reference-local, 1-based. Insertion coordinate is the reference base "
    "immediately before the insertion (0 means before the first base). Copied "
    "segment start and end are 1-based and inclusive. "
    "The copied segment is immediately before the insertion when it ends at the "
    "insertion coordinate, and immediately after when it starts at the following "
    "reference base. Before/after describes this coordinate representation, not "
    "a biological direction of copying."
)

FILTER_REASON_LABELS: dict[str, tuple[str, str]] = {
    "PARTIAL_OBSERVATION": (
        "PARTIAL",
        "The insertion is observed at a read edge, so the full ITD cannot be reconstructed.",
    ),
    "ONLY_CONFLICTING_MATE_EVIDENCE": (
        "CONFLICT",
        "Only conflicting mate evidence was observed for this candidate.",
    ),
    "ONLY_UNRESOLVED_EVIDENCE": (
        "UNRESOLVED",
        "Only unresolved evidence was observed for this candidate.",
    ),
    "AMBIGUOUS_EVIDENCE_DOMINATES": (
        "AMBIGUOUS",
        "Conflicting and unresolved evidence exceeds mutant and wild-type evidence.",
    ),
    "LOW_MUTANT_FRAGMENT_COUNT": (
        "LOW-SUPPORT",
        "Mutant-supporting fragment count is below the configured minimum.",
    ),
    "LOW_INFORMATIVE_FRAGMENT_COUNT": (
        "LOW-DEPTH",
        "Informative fragment count is below the configured minimum.",
    ),
    "LOW_MUTANT_FRAGMENT_FRACTION": (
        "LOW-FRACTION",
        "Observed mutant-fragment fraction is below the configured minimum.",
    ),
    "DIRECTION_BIAS": (
        "BIASED",
        "Directional mutant fractions are more imbalanced than the configured limit.",
    ),
}


class OptionalDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults only for options with concrete defaults."""

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help
        if help_text is None:
            help_text = ""
        if (
            "%(default)" not in help_text
            and action.default is not argparse.SUPPRESS
            and action.default is not None
            and not action.required
        ):
            help_text += " (default: %(default)s)"
        return help_text


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="itdiscover",
        description="Discover FLT3 ITDs from amplicon sequencing of AML samples.",
        formatter_class=OptionalDefaultsHelpFormatter,
        allow_abbrev=False,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference amplicon FASTA file containing exactly one sequence.",
    )
    parser.add_argument(
        "--r1",
        required=True,
        help="Forward-read FASTQ file.",
    )
    parser.add_argument(
        "--r2",
        required=True,
        help="Reverse-read FASTQ file.",
    )
    parser.add_argument(
        "--forward-primer",
        required=True,
        help="Forward primer sequence trimmed from the 5' end of R1 reads.",
    )
    parser.add_argument(
        "--reverse-primer",
        required=True,
        help=(
            "Reverse primer sequence as it occurs at the 5' end of raw R2; "
            "its reverse complement is trimmed from oriented R2 reads."
        ),
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=100,
        help="Minimum read length after trimming terminal Ns.",
    )
    parser.add_argument(
        "--min-mean-quality",
        type=float,
        default=30,
        help="Minimum mean Phred quality score per read.",
    )
    parser.add_argument(
        "--min-insert-length",
        type=_positive_int,
        default=6,
        help="Minimum insertion length to consider.",
    )
    parser.add_argument(
        "--min-copied-segment-length",
        type=_positive_int,
        help=(
            "Minimum copied reference-segment length; defaults to "
            "--min-insert-length."
        ),
    )
    parser.add_argument(
        "--require-in-frame",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require fully observed insertion length to be divisible by three. "
            "Use --no-require-in-frame to retain out-of-frame candidates."
        ),
    )
    parser.add_argument(
        "--max-copy-mismatch-rate",
        type=_fraction,
        default=0.0,
        help=(
            "Maximum copied-segment mismatches divided by copied-segment "
            "length; 0 is exact mode. ITD detection happens before any "
            "optional minor-allele consolidation."
        ),
    )
    consolidation_group = parser.add_argument_group(
        "advanced minor-allele consolidation"
    )
    consolidation_group.add_argument(
        "--consolidate-minor-itd-variants",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Consolidate weak, directly compatible ITD observations into a "
            "dominant allele. This is opt-in because it changes allele identity."
        ),
    )
    consolidation_group.add_argument(
        "--consolidation-max-allele-mismatch-rate",
        type=_fraction,
        default=0.125,
        help=(
            "Maximum complete-ALT positional mismatches divided by insertion "
            "length for an already-detected minor allele and its anchor. This "
            "does not affect whether an insertion qualifies as an ITD."
        ),
    )
    consolidation_group.add_argument(
        "--consolidation-max-breakpoint-shift-rate",
        type=_fraction,
        default=1.0,
        help=(
            "Maximum absolute breakpoint shift divided by insertion length "
            "between a minor allele and its anchor."
        ),
    )
    consolidation_group.add_argument(
        "--consolidation-max-minor-support-ratio",
        type=_fraction,
        default=0.05,
        help=(
            "Maximum evidence-passing raw fragment support of a minor allele "
            "relative to its anchor."
        ),
    )
    consolidation_group.add_argument(
        "--consolidation-min-anchor-fragment-count",
        type=_positive_int,
        default=3,
        help=(
            "Minimum evidence-passing raw fragment support required for an "
            "anchor allele."
        ),
    )
    parser.add_argument(
        "--min-mutant-fragment-count",
        type=_positive_int,
        default=3,
        help="Minimum mutant fragment count required to pass filtering.",
    )
    parser.add_argument(
        "--min-informative-fragment-count",
        type=_non_negative_int,
        default=10,
        help="Minimum mutant-plus-wild-type fragment count required to pass.",
    )
    parser.add_argument(
        "--min-mutant-fragment-fraction",
        type=_fraction,
        default=0.01,
        help=(
            "Minimum observed mutant/(mutant + wild type) fragment fraction "
            "required to pass; this is not a validated VAF or allelic ratio."
        ),
    )
    parser.add_argument(
        "--max-directional-mutant-fraction-share",
        type=_direction_fraction,
        default=0.90,
        help=(
            "Largest allowed share of the two direction-specific mutant "
            "fractions attributable to R1 or R2."
        ),
    )
    parser.add_argument(
        "--min-directional-opportunities",
        type=_positive_int,
        default=5,
        help=(
            "Minimum callable junction opportunities required in each "
            "direction before direction-bias filtering."
        ),
    )
    parser.add_argument(
        "--min-alignment-identity",
        type=_fraction,
        default=0.90,
        help="Minimum identity across read bases aligned to reference bases.",
    )
    parser.add_argument(
        "--min-on-target-fraction",
        type=_fraction,
        default=0.80,
        help="Minimum fraction of the shorter sequence aligned to the target.",
    )
    parser.add_argument(
        "--min-alignment-score",
        type=float,
        help="Optional minimum raw pairwise alignment score.",
    )
    parser.add_argument(
        "--reject-ambiguous-alignments",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Reject reads whose optimal alignments imply different normalized "
            "sequence events. Equivalent insertion-gap placements are retained."
        ),
    )
    parser.add_argument(
        "--min-junction-anchor-quality",
        type=_non_negative_int,
        default=30,
        help="Minimum Phred quality for every junction-anchor base.",
    )
    parser.add_argument(
        "--min-insert-mean-quality",
        type=_non_negative_float,
        default=30.0,
        help="Minimum mean Phred quality across inserted bases.",
    )
    parser.add_argument(
        "--min-insert-base-quality",
        type=_non_negative_int,
        default=15,
        help="Minimum Phred quality permitted for any inserted base.",
    )
    parser.add_argument(
        "--junction-flank-size",
        type=_positive_int,
        default=3,
        help="High-quality read bases required on each side of an insertion.",
    )
    parser.add_argument(
        "--sample-id",
        help="Sample identifier shown in reports; defaults to the R1 filename stem.",
    )
    parser.add_argument(
        "--min-usable-fragment-count",
        type=_positive_int,
        default=10,
        help="Minimum post-preprocessing fragments required for sample QC.",
    )
    parser.add_argument(
        "--min-qc-reads-per-direction",
        type=_non_negative_int,
        default=1,
        help="Minimum alignment-passing R1 and R2 reads required for sample QC.",
    )
    parser.add_argument(
        "--min-alignment-pass-fraction",
        type=_fraction,
        default=0.80,
        help="Minimum fraction of preprocessed reads passing alignment filters.",
    )
    parser.add_argument(
        "--min-median-interbase-coverage",
        type=_non_negative_int,
        default=10,
        help="Minimum median fragment coverage across target inter-base sites.",
    )
    parser.add_argument(
        "--min-primer-retention-fraction",
        type=_fraction,
        default=0.80,
        help="Minimum per-direction primer retention when a primer is configured.",
    )
    parser.add_argument(
        "--output",
        type=_html_output_path,
        help="Optional path for a concise HTML result summary.",
    )
    parser.add_argument(
        "--output-tsv",
        type=_tsv_output_path,
        help="Optional path for a TSV summary of called ITDs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the ITDiscover CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return _run_call_command(args)
    except PrimerOrientationError as error:
        _write_analysis_error_reports(args, error)
        print(f"{parser.prog}: error: {error}", file=sys.stderr)
        return 2
    except Exception as error:
        _write_analysis_error_reports(args, error)
        raise


def _run_call_command(args: argparse.Namespace) -> int:
    reference_id, reference = _read_single_sequence_fasta_record(
        Path(args.reference)
    )
    reference_sha256 = hashlib.sha256(reference.encode("ascii")).hexdigest()
    fragments = list(read_paired_fastq(args.r1, args.r2))
    trimming = _build_trim_settings(args)
    validate_primer_orientations(fragments, trimming)
    preprocessing_result = preprocess_fragments_with_metrics(
        fragments,
        min_length=args.min_read_length,
        min_mean_quality=args.min_mean_quality,
        trimming=trimming,
    )
    processed_reads = list(preprocessing_result.reads)
    alignment_filters = AlignmentEvidenceFilter(
        min_identity=args.min_alignment_identity,
        min_on_target_fraction=args.min_on_target_fraction,
        min_score=args.min_alignment_score,
        reject_ambiguous=args.reject_ambiguous_alignments,
    )
    unfiltered_alignments = [
        align_read_to_reference(
            read,
            reference,
            detect_ambiguous_events=alignment_filters.reject_ambiguous,
        )
        for read in processed_reads
    ]
    alignments = [
        alignment
        for alignment in unfiltered_alignments
        if passes_alignment_evidence_filters(alignment, alignment_filters)
    ]
    sample_alignment_metrics = alignment_metrics(len(processed_reads), alignments)
    sample_coverage_metrics = coverage_metrics(alignments, len(reference))
    qc_thresholds = SampleQCThresholds(
        min_usable_fragment_count=args.min_usable_fragment_count,
        min_passing_reads_per_direction=args.min_qc_reads_per_direction,
        min_alignment_pass_fraction=args.min_alignment_pass_fraction,
        min_median_interbase_coverage=args.min_median_interbase_coverage,
        min_primer_retention_fraction=args.min_primer_retention_fraction,
    )
    insertion_filters = InsertionEvidenceFilter(
        min_junction_anchor_quality=args.min_junction_anchor_quality,
        min_insert_mean_quality=args.min_insert_mean_quality,
        min_insert_base_quality=args.min_insert_base_quality,
        junction_flank_size=args.junction_flank_size,
    )
    filters = ITDFilter(
        min_mutant_fragment_count=args.min_mutant_fragment_count,
        min_informative_fragment_count=args.min_informative_fragment_count,
        min_observed_mutant_fragment_fraction=(
            args.min_mutant_fragment_fraction
        ),
        max_directional_mutant_fraction_share=(
            args.max_directional_mutant_fraction_share
        ),
        min_directional_opportunities=args.min_directional_opportunities,
    )
    consolidation = _build_consolidation_settings(args)
    min_copied_segment_length = (
        args.min_insert_length
        if args.min_copied_segment_length is None
        else args.min_copied_segment_length
    )
    if args.max_copy_mismatch_rate == 0:
        calls, representatives = call_exact_itds_with_representatives(
            alignments,
            reference,
            min_insert_length=args.min_insert_length,
            min_copied_segment_length=min_copied_segment_length,
            require_in_frame=args.require_in_frame,
            filters=filters,
            evidence_filter=insertion_filters,
            consolidation=consolidation,
        )
    else:
        calls, representatives = call_fuzzy_itds_with_representatives(
            alignments,
            reference,
            max_copy_mismatch_rate=args.max_copy_mismatch_rate,
            min_insert_length=args.min_insert_length,
            min_copied_segment_length=min_copied_segment_length,
            require_in_frame=args.require_in_frame,
            filters=filters,
            evidence_filter=insertion_filters,
            consolidation=consolidation,
        )
    sample_result = build_sample_result(
        sample_id=_sample_id(args),
        calls=calls,
        preprocessing=preprocessing_result.metrics,
        alignment=sample_alignment_metrics,
        coverage=sample_coverage_metrics,
        thresholds=qc_thresholds,
    )
    if args.output:
        _write_html_report(
            args.output,
            calls,
            representatives,
            sample_result=sample_result,
            show_sample_id=bool(args.sample_id),
        )
    if args.output_tsv:
        _write_tsv_call_report(
            args.output_tsv,
            calls,
            max_copy_mismatch_rate=args.max_copy_mismatch_rate,
            min_mutant_fragment_count=filters.min_mutant_fragment_count,
            min_informative_fragment_count=filters.min_informative_fragment_count,
            min_mutant_fragment_fraction=(
                filters.min_observed_mutant_fragment_fraction
            ),
            max_directional_mutant_fraction_share=(
                filters.max_directional_mutant_fraction_share
            ),
            min_directional_opportunities=filters.min_directional_opportunities,
            alignment_filters=alignment_filters,
            insertion_filters=insertion_filters,
            sample_result=sample_result,
            qc_thresholds=qc_thresholds,
            min_insert_length=args.min_insert_length,
            min_copied_segment_length=min_copied_segment_length,
            require_in_frame=args.require_in_frame,
            consolidation=consolidation,
            reference_id=reference_id,
            reference_length=len(reference),
            reference_sha256=reference_sha256,
        )
    return 0


def _html_output_path(value: str) -> Path:
    path = Path(value)
    if path.suffix.lower() != ".html":
        raise argparse.ArgumentTypeError("output path must end with .html")
    return path


def _tsv_output_path(value: str) -> Path:
    path = Path(value)
    if path.suffix.lower() != ".tsv":
        raise argparse.ArgumentTypeError("output path must end with .tsv")
    return path


def _non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must not be negative")
    return parsed


def _non_negative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must not be negative")
    return parsed


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def _fraction(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError("value must be between 0 and 1")
    return parsed


def _direction_fraction(value: str) -> float:
    parsed = _fraction(value)
    if parsed < 0.5:
        raise argparse.ArgumentTypeError("value must be between 0.5 and 1")
    return parsed


def _build_trim_settings(args: argparse.Namespace) -> ReadTrimSettings:
    return ReadTrimSettings(
        forward_primer=args.forward_primer,
        reverse_primer=args.reverse_primer,
    )


def _build_consolidation_settings(
    args: argparse.Namespace,
) -> ITDConsolidationSettings:
    return ITDConsolidationSettings(
        enabled=args.consolidate_minor_itd_variants,
        max_allele_mismatch_rate=(
            args.consolidation_max_allele_mismatch_rate
        ),
        max_breakpoint_shift_rate=(
            args.consolidation_max_breakpoint_shift_rate
        ),
        max_minor_to_anchor_support_ratio=(
            args.consolidation_max_minor_support_ratio
        ),
        min_anchor_fragment_count=(
            args.consolidation_min_anchor_fragment_count
        ),
    )


def _sample_id(args: argparse.Namespace) -> str:
    if args.sample_id:
        return args.sample_id
    name = Path(args.r1).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for suffix in ("_R1", "-R1", ".R1"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def _write_analysis_error_reports(
    args: argparse.Namespace,
    error: Exception,
) -> None:
    result = error_sample_result(_sample_id(args), error)
    reference_id, reference_length, reference_sha256 = (
        _reference_report_metadata(Path(args.reference))
    )
    qc_thresholds = SampleQCThresholds(
        min_usable_fragment_count=args.min_usable_fragment_count,
        min_passing_reads_per_direction=args.min_qc_reads_per_direction,
        min_alignment_pass_fraction=args.min_alignment_pass_fraction,
        min_median_interbase_coverage=args.min_median_interbase_coverage,
        min_primer_retention_fraction=args.min_primer_retention_fraction,
    )
    consolidation = _build_consolidation_settings(args)
    if args.output:
        try:
            _write_html_report(
                args.output,
                [],
                [],
                sample_result=result,
                show_sample_id=bool(args.sample_id),
            )
        except Exception:
            pass
    if args.output_tsv:
        try:
            _write_tsv_call_report(
                args.output_tsv,
                [],
                max_copy_mismatch_rate=args.max_copy_mismatch_rate,
                min_mutant_fragment_count=args.min_mutant_fragment_count,
                min_informative_fragment_count=args.min_informative_fragment_count,
                min_mutant_fragment_fraction=(
                    args.min_mutant_fragment_fraction
                ),
                sample_result=result,
                qc_thresholds=qc_thresholds,
                min_insert_length=args.min_insert_length,
                min_copied_segment_length=(
                    args.min_insert_length
                    if args.min_copied_segment_length is None
                    else args.min_copied_segment_length
                ),
                require_in_frame=args.require_in_frame,
                consolidation=consolidation,
                reference_id=reference_id,
                reference_length=reference_length,
                reference_sha256=reference_sha256,
            )
        except Exception:
            pass


def _format_filter_reasons(call: ITDCall) -> str:
    return "." if not call.filter_reasons else ";".join(call.filter_reasons)


def _format_consolidated_members(call: ITDCall) -> str:
    if not call.consolidated_members:
        return "."
    return " | ".join(
        (
            f"start={member.allele.start},"
            f"sequence={member.allele.sequence},"
            f"fragments={member.fragment_count},"
            f"allele_mismatches={member.allele_mismatches},"
            f"allele_mismatch_rate={member.allele_mismatch_rate:.6f},"
            f"breakpoint_shift={member.breakpoint_shift},"
            f"breakpoint_shift_rate={member.breakpoint_shift_rate:.6f},"
            f"reason={member.reason}"
        )
        for member in call.consolidated_members
    )


def _format_directional_evidence(
    mutant_count: int,
    opportunity_count: int | None,
) -> str:
    """Format one direction's mutant count, opportunities, and fraction."""
    if opportunity_count in (None, 0):
        return "not evaluable (0 opportunities)"
    return (
        f"{mutant_count / opportunity_count:.1%} "
        f"({mutant_count}/{opportunity_count} opportunities)"
    )


def _read_single_sequence_fasta(path: Path) -> str:
    """Return the only sequence in a FASTA file."""
    _, sequence = _read_single_sequence_fasta_record(path)
    return sequence


def _reference_report_metadata(path: Path) -> tuple[str, int | None, str | None]:
    try:
        reference_id, sequence = _read_single_sequence_fasta_record(path)
    except Exception:
        return path.name, None, None
    return (
        reference_id,
        len(sequence),
        hashlib.sha256(sequence.encode("ascii")).hexdigest(),
    )


def _read_single_sequence_fasta_record(path: Path) -> tuple[str, str]:
    with path.open(mode="rt", encoding="utf-8") as handle:
        records = list(_iter_fasta_records(handle))
    if not records:
        raise ValueError("reference FASTA does not contain a sequence")
    if len(records) > 1:
        raise ValueError("reference FASTA must contain exactly one sequence")
    return records[0]


def _iter_fasta_records(handle: TextIO) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_name: str | None = None
    current_parts: list[str] = []

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                records.append((current_name, "".join(current_parts)))
            current_name = line[1:].strip()
            if not current_name:
                raise ValueError("reference FASTA header must not be empty")
            current_parts = []
            continue
        if current_name is None:
            raise ValueError("reference FASTA sequence must follow a header")
        current_parts.append(line)

    if current_name is not None:
        records.append((current_name, "".join(current_parts)))
    return records


def _iter_fasta_sequences(handle: TextIO) -> list[str]:
    """Return FASTA sequences without their record identifiers."""
    return [sequence for _, sequence in _iter_fasta_records(handle)]


def _write_html_report(
    path: Path,
    calls: list[ITDCall],
    representatives: list[UniqueSupportRepresentative],
    *,
    sample_result: SampleResult | None = None,
    show_sample_id: bool = False,
) -> None:
    representatives_by_allele: dict[
        object, list[UniqueSupportRepresentative]
    ] = {}
    for representative in representatives:
        representatives_by_allele.setdefault(
            representative.canonical_allele,
            [],
        ).append(representative)

    show_calls = sample_result is None or sample_result.outcome == "ITD detected"
    ordered_calls = sorted(
        (call for call in calls if show_calls and call.status == "PASS"),
        key=_html_call_sort_key,
    )
    sections = [
        _render_html_call_section(
            call,
            _best_representative(
                representatives_by_allele.get(call.canonical_allele, [])
            ),
        )
        for call in ordered_calls
    ]
    filtered_calls = sorted(
        (call for call in calls if call.status != "PASS"),
        key=_html_call_sort_key,
    )
    filtered_variants = _render_html_filtered_variants(
        filtered_calls,
        representatives_by_allele,
    )

    sample_summary = _render_html_sample_summary(sample_result, show_sample_id)
    empty_state = "" if sections else _render_html_empty_state(sample_result)

    document = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>ITDiscover Report</title>
  <style>
    :root {
      color-scheme: light;
      --ink: #17212b;
      --muted: #607080;
      --line: #dce3e8;
      --panel: #f5f8fa;
      --pass: #13734b;
      --pass-bg: #e8f5ee;
      --warn: #8a5b00;
      --warn-bg: #fff5d6;
      --fail: #a33131;
      --fail-bg: #fdecec;
    }
    body {
      margin: 0;
      padding: 40px 20px;
      color: var(--ink);
      font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
      line-height: 1.45;
      background: #eef2f5;
    }
    main {
      max-width: 860px;
      margin: 0 auto;
    }
    h2, h3 { margin: 0; }
    .sample-result, .itd, .empty-state {
      background: white;
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 22px;
      margin-bottom: 18px;
      box-shadow: 0 1px 2px rgb(23 33 43 / 4%);
    }
    .sample-name {
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.07em;
    }
    .outcome {
      margin: 4px 0 16px;
      font-size: 26px;
      line-height: 1.2;
    }
    .status-row {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-bottom: 16px;
    }
    .status {
      display: inline-flex;
      align-items: center;
      gap: 5px;
      font-size: 13px;
      font-weight: 650;
    }
    .status-value, .filter-reason {
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 3px 9px;
      border: 1px solid transparent;
      font-size: 12px;
      font-weight: 750;
      letter-spacing: 0.02em;
      line-height: 1.35;
    }
    .status-value--success {
      color: #166534;
      background: #bbf7d0;
      border-color: #4ade80;
    }
    .status-value--failure {
      color: #991b1b;
      background: #fecaca;
      border-color: #f87171;
    }
    .filter-reasons {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
    }
    .filter-reason {
      color: #991b1b;
      background: #fecaca;
      border-color: #f87171;
      padding: 1px 6px;
      font-size: 10px;
      letter-spacing: 0.04em;
      cursor: help;
    }
    .summary {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(175px, 1fr));
      gap: 14px 22px;
      margin: 16px 0 0;
    }
    .summary dt {
      margin: 0 0 4px;
      color: var(--muted);
      font-size: 12px;
    }
    .summary dd {
      margin: 0;
      font-size: 15px;
      font-weight: 600;
      overflow-wrap: anywhere;
    }
    .sequence {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 17px !important;
      letter-spacing: 0.03em;
    }
    .note, .alert {
      margin: 16px 0 0;
      font-size: 14px;
    }
    .note { color: var(--muted); }
    .alert {
      background: var(--panel);
      border-radius: 7px;
      padding: 10px 12px;
    }
    .itd .summary { margin-top: 0; }
    .alignment {
      margin-top: 18px;
      border-top: 1px solid var(--line);
      padding-top: 14px;
    }
    .alignment-block {
      margin-top: 10px;
      overflow-x: auto;
      border: 1px solid var(--line);
      border-radius: 7px;
      padding: 10px;
      padding-top: 32px;
      background: var(--panel);
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      position: relative;
    }
    .alignment-row { white-space: pre; }
    .alignment-label { color: var(--muted); }
    .alignment-ruler {
      position: absolute;
      top: 8px;
      left: 10px;
      height: 16px;
      color: var(--muted);
      white-space: nowrap;
    }
    .position-marker {
      position: absolute;
      top: 0;
      color: var(--ink);
      font-weight: 700;
      transform: translateX(-50%);
    }
    .alignment-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px 14px;
      margin-top: 10px;
      color: var(--muted);
      font-size: 12px;
    }
    .legend-item { display: inline-flex; align-items: center; gap: 5px; }
    .legend-chip {
      border-radius: 3px;
      padding: 1px 5px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 11px;
      font-weight: 700;
    }
    .tandem-region {
      background: #dcfce7;
      color: #166534;
      font-weight: 700;
    }
    .inserted-region {
      background: #dbeafe;
      color: #12315d;
      font-weight: 700;
    }
    .spacer-region {
      background: #fef3c7;
      color: #92400e;
      font-weight: 700;
    }
    .diff, .insert-mismatch {
      background: #fee2e2;
      color: #b91c1c;
      font-weight: 700;
    }
    .filtered-variants {
      margin-bottom: 18px;
      padding: 18px 22px;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: white;
      box-shadow: 0 1px 2px rgb(23 33 43 / 4%);
    }
    .filtered-variants summary {
      cursor: pointer;
      color: var(--muted);
      font-weight: 700;
    }
    .filtered-variants-content {
      display: grid;
      gap: 14px;
      margin-top: 14px;
      min-width: 0;
    }
    .filtered-variant {
      margin: 0;
      min-width: 0;
    }
    @media (max-width: 520px) {
      body { padding: 24px 12px; }
      .sample-result, .itd, .empty-state { padding: 18px; }
    }
  </style>
</head>
<body>
  <main>
    __SAMPLE_SUMMARY__
    __EMPTY_STATE__
    __SECTIONS__
    __FILTERED_VARIANTS__
  </main>
</body>
</html>
"""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        document.replace("__SAMPLE_SUMMARY__", sample_summary)
        .replace("__EMPTY_STATE__", empty_state)
        .replace("__SECTIONS__", "\n".join(sections))
        .replace("__FILTERED_VARIANTS__", filtered_variants),
        encoding="utf-8",
    )


def _write_tsv_call_report(
    path: Path,
    calls: list[ITDCall],
    *,
    max_copy_mismatch_rate: float,
    min_mutant_fragment_count: int,
    min_informative_fragment_count: int,
    min_mutant_fragment_fraction: float,
    max_directional_mutant_fraction_share: float = 0.90,
    min_directional_opportunities: int = 5,
    alignment_filters: AlignmentEvidenceFilter | None = None,
    insertion_filters: InsertionEvidenceFilter | None = None,
    sample_result: SampleResult | None = None,
    qc_thresholds: SampleQCThresholds | None = None,
    min_insert_length: int = 6,
    min_copied_segment_length: int = 6,
    require_in_frame: bool = True,
    consolidation: ITDConsolidationSettings = ITDConsolidationSettings(),
    reference_id: str | None = None,
    reference_length: int | None = None,
    reference_sha256: str | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open(mode="wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Call",
                "Status",
                "Filter Reasons",
                "Mode",
                "Max Copy Mismatch Rate",
                "Insertion After Reference Base (1-based; 0=before first)",
                "Copied Segment Start (1-based)",
                "Copied Segment End (1-based, inclusive)",
                "Copied Segment Sequence",
                "Spacer Prefix",
                "Spacer Suffix",
                "Insertion Sequence",
                "Read-Edge Observation",
                "Mutant Fragment Count",
                "R1 Mutant Count",
                "R1 Opportunity Count",
                "R1 Mutant Fraction",
                "R2 Mutant Count",
                "R2 Opportunity Count",
                "R2 Mutant Fraction",
                "Informative Fragment Count",
                "Observed Mutant-fragment Fraction",
                "Mutant/Informative Fragments",
                "Min Mutant Fragment Count",
                "Min Informative Fragment Count",
                "Min Mutant-fragment Fraction",
                "Max Directional Mutant-fraction Share",
                "Min Opportunities per Direction",
                "Min Alignment Identity",
                "Min On-target Fraction",
                "Min Alignment Score",
                "Reject Ambiguous Alignments",
                "Min Junction Anchor Quality",
                "Min Insert Mean Quality",
                "Min Insert Base Quality",
                "Junction Flank Size",
                "Sample ID",
                "Analysis Status",
                "QC Status",
                "Outcome",
                "QC Reasons",
                "Analysis Error",
                "Input Fragment Count",
                "Input Read Count",
                "Forward Primer-retained Reads",
                "Reverse Primer-retained Reads",
                "Primer-failed Read Count",
                "Length-failed Read Count",
                "Quality-failed Read Count",
                "Preprocessing-passing Read Count",
                "Preprocessing-passing Forward Reads",
                "Preprocessing-passing Reverse Reads",
                "Usable Fragment Count",
                "Alignment-attempted Read Count",
                "Alignment-passing Read Count",
                "Alignment-passing Forward Reads",
                "Alignment-passing Reverse Reads",
                "Alignment-passing Fragment Count",
                "Alignment Pass Fraction",
                "Minimum Inter-base Coverage",
                "Median Inter-base Coverage",
                "Maximum Inter-base Coverage",
                "Passing Call Count",
                "Filtered Candidate Count",
                "QC Min Usable Fragments",
                "QC Min Reads per Direction",
                "QC Min Alignment Pass Fraction",
                "QC Min Median Inter-base Coverage",
                "QC Min Primer Retention Fraction",
                "Concordant Fragment Count",
                "Single-mate Fragment Count",
                "Conflicting Fragment Count",
                "Unresolved Fragment Count",
                "Wild-type Fragment Count",
                "Not-informative Fragment Count",
                "Reference FASTA Header",
                "Reference Length",
                "Reference Sequence SHA-256",
                "Coordinate Convention",
                "Copied Segment Location",
                "Min Insert Length",
                "Min Copied Segment Length",
                "Require In-frame Insertions",
                "Minor-variant Consolidation Enabled",
                "Consolidation Max Allele Mismatch Rate",
                "Consolidation Max Breakpoint Shift Rate",
                "Consolidation Max Minor/Anchor Support Ratio",
                "Consolidation Min Anchor Fragment Count",
                "Consolidated Minor Allele Count",
                "Consolidated Minor Raw Fragment Support",
                "Consolidated Minor Alleles",
            ]
        )
        mode = "exact" if max_copy_mismatch_rate == 0 else "fuzzy"
        for index, call in enumerate(calls, start=1):
            writer.writerow(
                [
                    index,
                    call.status,
                    _format_filter_reasons(call),
                    mode,
                    f"{max_copy_mismatch_rate:.6f}",
                    call.itd.insertion.start + 1,
                    call.itd.copied_segment_start + 1,
                    call.itd.copied_segment_end + 1,
                    call.itd.copied_segment_sequence,
                    call.itd.spacer_prefix or "-",
                    call.itd.spacer_suffix or "-",
                    call.itd.insertion.sequence,
                    "Yes" if call.itd.is_partial_observation else "No",
                    call.mutant_fragment_count,
                    call.r1_mutant_count,
                    call.r1_opportunity_count,
                    (
                        f"{call.r1_mutant_fraction:.6f}"
                        if call.r1_mutant_fraction is not None
                        else "."
                    ),
                    call.r2_mutant_count,
                    call.r2_opportunity_count,
                    (
                        f"{call.r2_mutant_fraction:.6f}"
                        if call.r2_mutant_fraction is not None
                        else "."
                    ),
                    call.informative_fragment_count,
                    f"{call.observed_mutant_fragment_fraction:.6f}",
                    (
                        f"{call.mutant_fragment_count}/"
                        f"{call.informative_fragment_count} informative fragments"
                    ),
                    min_mutant_fragment_count,
                    min_informative_fragment_count,
                    f"{min_mutant_fragment_fraction:.6f}",
                    f"{max_directional_mutant_fraction_share:.6f}",
                    min_directional_opportunities,
                    (
                        f"{alignment_filters.min_identity:.6f}"
                        if alignment_filters is not None
                        else "."
                    ),
                    (
                        f"{alignment_filters.min_on_target_fraction:.6f}"
                        if alignment_filters is not None
                        else "."
                    ),
                    (
                        alignment_filters.min_score
                        if alignment_filters is not None
                        and alignment_filters.min_score is not None
                        else "."
                    ),
                    (
                        "Yes"
                        if alignment_filters is not None
                        and alignment_filters.reject_ambiguous
                        else "No"
                    ),
                    (
                        insertion_filters.min_junction_anchor_quality
                        if insertion_filters is not None
                        else "."
                    ),
                    (
                        f"{insertion_filters.min_insert_mean_quality:.6f}"
                        if insertion_filters is not None
                        else "."
                    ),
                    (
                        insertion_filters.min_insert_base_quality
                        if insertion_filters is not None
                        else "."
                    ),
                    (
                        insertion_filters.junction_flank_size
                        if insertion_filters is not None
                        else "."
                    ),
                    *_sample_tsv_values(sample_result, qc_thresholds),
                    call.concordant_fragment_count,
                    call.single_mate_fragment_count,
                    call.conflicting_fragment_count,
                    call.unresolved_fragment_count,
                    call.wild_type_fragment_count,
                    call.not_informative_fragment_count,
                    reference_id or ".",
                    reference_length if reference_length is not None else ".",
                    reference_sha256 or ".",
                    COORDINATE_CONVENTION,
                    f"immediately {call.itd.copied_segment_location} insertion",
                    min_insert_length,
                    min_copied_segment_length,
                    "Yes" if require_in_frame else "No",
                    "Yes" if consolidation.enabled else "No",
                    f"{consolidation.max_allele_mismatch_rate:.6f}",
                    f"{consolidation.max_breakpoint_shift_rate:.6f}",
                    f"{consolidation.max_minor_to_anchor_support_ratio:.6f}",
                    consolidation.min_anchor_fragment_count,
                    len(call.consolidated_members),
                    call.consolidated_minor_fragment_count,
                    _format_consolidated_members(call),
                ]
            )
        if not calls:
            empty_call_values: list[object] = ["."] * 36
            empty_call_values[3] = mode
            empty_call_values[4] = f"{max_copy_mismatch_rate:.6f}"
            empty_call_values[23] = min_mutant_fragment_count
            empty_call_values[24] = min_informative_fragment_count
            empty_call_values[25] = f"{min_mutant_fragment_fraction:.6f}"
            empty_call_values[26] = f"{max_directional_mutant_fraction_share:.6f}"
            empty_call_values[27] = min_directional_opportunities
            if alignment_filters is not None:
                empty_call_values[28] = f"{alignment_filters.min_identity:.6f}"
                empty_call_values[29] = (
                    f"{alignment_filters.min_on_target_fraction:.6f}"
                )
                empty_call_values[30] = (
                    alignment_filters.min_score
                    if alignment_filters.min_score is not None
                    else "."
                )
                empty_call_values[31] = (
                    "Yes" if alignment_filters.reject_ambiguous else "No"
                )
            if insertion_filters is not None:
                empty_call_values[32] = (
                    insertion_filters.min_junction_anchor_quality
                )
                empty_call_values[33] = (
                    f"{insertion_filters.min_insert_mean_quality:.6f}"
                )
                empty_call_values[34] = insertion_filters.min_insert_base_quality
                empty_call_values[35] = insertion_filters.junction_flank_size
            writer.writerow(
                empty_call_values
                + _sample_tsv_values(sample_result, qc_thresholds)
                + [".", ".", ".", ".", ".", "."]
                + [
                    reference_id or ".",
                    reference_length if reference_length is not None else ".",
                    reference_sha256 or ".",
                    COORDINATE_CONVENTION,
                    ".",
                ]
                + [
                    min_insert_length,
                    min_copied_segment_length,
                    "Yes" if require_in_frame else "No",
                    "Yes" if consolidation.enabled else "No",
                    f"{consolidation.max_allele_mismatch_rate:.6f}",
                    f"{consolidation.max_breakpoint_shift_rate:.6f}",
                    f"{consolidation.max_minor_to_anchor_support_ratio:.6f}",
                    consolidation.min_anchor_fragment_count,
                    ".",
                    ".",
                    ".",
                ]
            )


def _sample_tsv_values(
    result: SampleResult | None,
    thresholds: SampleQCThresholds | None,
) -> list[object]:
    if result is None:
        result_values: list[object] = ["."] * 28
    else:
        preprocessing = result.preprocessing
        alignment = result.alignment
        coverage = result.coverage
        result_values = [
            result.sample_id,
            result.analysis_status,
            result.qc_status,
            result.outcome,
            ";".join(result.qc_reasons) or ".",
            result.error_message or ".",
            preprocessing.input_fragment_count if preprocessing else ".",
            preprocessing.input_read_count if preprocessing else ".",
            (
                preprocessing.primer_retained_forward_reads
                if preprocessing
                and preprocessing.primer_retained_forward_reads is not None
                else "."
            ),
            (
                preprocessing.primer_retained_reverse_reads
                if preprocessing
                and preprocessing.primer_retained_reverse_reads is not None
                else "."
            ),
            preprocessing.primer_failed_read_count if preprocessing else ".",
            preprocessing.length_failed_read_count if preprocessing else ".",
            preprocessing.quality_failed_read_count if preprocessing else ".",
            preprocessing.passing_read_count if preprocessing else ".",
            preprocessing.passing_forward_read_count if preprocessing else ".",
            preprocessing.passing_reverse_read_count if preprocessing else ".",
            preprocessing.usable_fragment_count if preprocessing else ".",
            alignment.attempted_read_count if alignment else ".",
            alignment.passing_read_count if alignment else ".",
            alignment.passing_forward_read_count if alignment else ".",
            alignment.passing_reverse_read_count if alignment else ".",
            alignment.passing_fragment_count if alignment else ".",
            f"{alignment.pass_fraction:.6f}" if alignment else ".",
            coverage.minimum if coverage else ".",
            f"{coverage.median:.1f}" if coverage else ".",
            coverage.maximum if coverage else ".",
            result.passing_call_count,
            result.filtered_candidate_count,
        ]
    threshold_values: list[object] = (
        ["."] * 5
        if thresholds is None
        else [
            thresholds.min_usable_fragment_count,
            thresholds.min_passing_reads_per_direction,
            f"{thresholds.min_alignment_pass_fraction:.6f}",
            thresholds.min_median_interbase_coverage,
            f"{thresholds.min_primer_retention_fraction:.6f}",
        ]
    )
    return result_values + threshold_values


def _render_html_empty_state(result: SampleResult | None) -> str:
    if result is not None and result.analysis_status == "error":
        message = "ITD calling did not complete; no biological outcome is available."
    elif result is not None and result.outcome == "indeterminate":
        message = (
            "No passing ITD result is available because sample QC did not pass."
        )
    else:
        message = "No passing ITD was detected."
    return f'<section class="empty-state">{html.escape(message)}</section>'


def _render_html_sample_summary(
    result: SampleResult | None,
    show_sample_id: bool,
) -> str:
    if result is None:
        return ""

    values: list[tuple[str, str]] = []
    if result.preprocessing is not None:
        values.append(
            ("Usable fragments", str(result.preprocessing.usable_fragment_count))
        )
    if result.alignment is not None:
        values.append(
            ("Alignment pass rate", f"{result.alignment.pass_fraction:.1%}")
        )
    if result.coverage is not None:
        values.append(
            ("Median coverage", f"{result.coverage.median:.1f}×")
        )

    value_html = "".join(
        f"<div><dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd></div>"
        for label, value in values
    )
    qc_reasons = (
        f'<p class="alert"><strong>QC:</strong> '
        f'{html.escape("; ".join(result.qc_reasons))}</p>'
        if result.qc_reasons
        else ""
    )
    error_message = (
        f'<p class="alert"><strong>Analysis error:</strong> '
        f'{html.escape(result.error_message)}</p>'
        if result.error_message
        else ""
    )
    call_word = "call" if result.passing_call_count == 1 else "calls"
    candidate_word = (
        "candidate" if result.filtered_candidate_count == 1 else "candidates"
    )
    analysis_status_class = (
        "success" if result.analysis_status == "complete" else "failure"
    )
    qc_status_class = "success" if result.qc_status == "pass" else "failure"
    sample_title_html = (
        f'<h2 class="outcome">{html.escape(result.sample_id)}</h2>'
        if show_sample_id
        else ""
    )
    return (
        '<section class="sample-result">'
        '<div class="sample-name">ITDiscover Report</div>'
        f"{sample_title_html}"
        '<div class="status-row">'
        '<span class="status"><span>Analysis Status:</span>'
        f'<span class="status-value status-value--{analysis_status_class}">'
        f'{html.escape(result.analysis_status.upper())}</span></span>'
        '<span class="status"><span>QC Status:</span>'
        f'<span class="status-value status-value--{qc_status_class}">'
        f'{html.escape(result.qc_status.upper())}</span></span>'
        "</div>"
        f'<dl class="summary">{value_html}</dl>'
        f"{qc_reasons}{error_message}"
        f'<p class="note">{result.passing_call_count} passing {call_word}; '
        f'{result.filtered_candidate_count} filtered {candidate_word}. '
        "Reference coordinates are 1-based. See the TSV output for full QC "
        "metrics, thresholds, and audit details.</p>"
        "</section>"
    )


def _html_call_sort_key(call: ITDCall) -> tuple[object, ...]:
    return (
        -call.mutant_fragment_count,
        call.itd.insertion.start,
        call.itd.copied_segment_start,
        call.itd.copied_segment_sequence,
        call.itd.spacer_prefix,
        call.itd.spacer_suffix,
        call.itd.insertion.sequence,
    )


def _render_html_filtered_variants(
    calls: list[ITDCall],
    representatives_by_allele: dict[object, list[UniqueSupportRepresentative]],
) -> str:
    if not calls:
        return ""
    variant_word = "variant" if len(calls) == 1 else "variants"
    cards = "".join(
        _render_html_call_section(
            call,
            _best_representative(
                representatives_by_allele.get(call.canonical_allele, [])
            ),
            filtered=True,
        )
        for call in calls
    )
    return (
        '<details class="filtered-variants">'
        f"<summary>Filtered {variant_word} ({len(calls)})</summary>"
        f'<div class="filtered-variants-content">{cards}</div>'
        "</details>"
    )


def _render_html_call_section(
    call: ITDCall,
    representative: UniqueSupportRepresentative | None,
    *,
    filtered: bool = False,
) -> str:
    insertion_label = "Insertion after reference base"
    insertion_position = str(call.itd.insertion.start + 1)
    if call.itd.insertion.start == -1:
        insertion_label = "Insertion position"
        insertion_position = "Before reference base 1"

    summary: tuple[tuple[str, str], ...] = (
        ("Inserted sequence", call.itd.insertion.sequence),
        (insertion_label, insertion_position),
        (
            "Copied reference bases",
            f"{call.itd.copied_segment_start + 1}–{call.itd.copied_segment_end + 1}",
        ),
        ("Copied sequence", call.itd.copied_segment_sequence),
        (
            'Observed mutant-fragment fraction',
            f"{call.observed_mutant_fragment_fraction:.1%}",
        ),
        ('Mutant fragments', str(call.mutant_fragment_count)),
        ('Informative fragments', str(call.informative_fragment_count)),
    )
    if filtered:
        summary += (("Filter reasons", ""),)
    summary_parts: list[str] = []
    for label, value in summary:
        css_class = (
            ' class="sequence"'
            if label in {"Inserted sequence", "Copied sequence"}
            else ""
        )
        if label == "Filter reasons":
            summary_parts.append(
                f"<div><dt>{html.escape(label)}</dt>"
                f'<dd class="filter-reasons">'
                f"{_render_html_filter_reason_labels(call)}</dd></div>"
            )
        else:
            summary_parts.append(
                f"<div><dt>{html.escape(label)}</dt>"
                f"<dd{css_class}>{html.escape(value)}</dd></div>"
            )
    summary_html = "".join(summary_parts)
    alignment_html = (
        _render_html_representative_alignment(representative)
        if representative is not None
        else ""
    )
    return (
        f'<section class="itd{" filtered-variant" if filtered else ""}">'
        f'<dl class="summary">{summary_html}</dl>'
        f"{alignment_html}"
        "</section>"
    )


def _render_html_filter_reason_labels(call: ITDCall) -> str:
    labels: list[str] = []
    for reason in call.filter_reasons:
        label, description = FILTER_REASON_LABELS.get(
            reason,
            (reason.replace("_", "-"), reason.replace("_", " ").capitalize()),
        )
        labels.append(
            f'<span class="filter-reason" title="'
            f'{html.escape(description, quote=True)}">{html.escape(label)}</span>'
        )
    return "".join(labels)


def _best_representative(
    representatives: list[UniqueSupportRepresentative],
) -> UniqueSupportRepresentative | None:
    if not representatives:
        return None
    return min(
        representatives,
        key=lambda representative: (
            representative.exact_support_count == 0,
            -representative.exact_support_count,
            -representative.support_count,
            representative.mismatches,
            representative.signature,
            representative.alignment.read_id,
        ),
    )


def _render_html_representative_alignment(
    representative: UniqueSupportRepresentative,
) -> str:
    alignment = representative.alignment
    reference_html = _highlight_alignment(
        alignment.aligned_reference,
        _reference_tandem_classes(alignment.aligned_reference, representative.itd),
    )
    read_html = _highlight_alignment(
        alignment.aligned_read,
        _alignment_difference_classes(alignment, representative.itd),
    )
    position_markers = _render_reference_position_markers(
        alignment.aligned_reference
    )
    return (
        '<section class="alignment">'
        '<div class="alignment-legend">'
        '<span class="legend-item"><span class="legend-chip tandem-region">C</span>copied reference</span>'
        '<span class="legend-item"><span class="legend-chip inserted-region">I</span>inserted</span>'
        '<span class="legend-item"><span class="legend-chip spacer-region">S</span>spacer</span>'
        '<span class="legend-item"><span class="legend-chip diff">M</span>mismatch</span>'
        '</div>'
        '<div class="alignment-block">'
        f'<div class="alignment-ruler"><span class="alignment-label">position   </span>{position_markers}</div>'
        '<div class="alignment-row"><span class="alignment-label">reference  </span>'
        f"{reference_html}</div>"
        '<div class="alignment-row"><span class="alignment-label">read       </span>'
        f"{read_html}</div>"
        "</div>"
        "</section>"
    )


def _highlight_alignment(sequence: str, classes: list[str | None]) -> str:
    fragments: list[str] = []
    for base, css_class in zip(sequence, classes, strict=True):
        escaped_base = html.escape(base)
        if css_class is not None:
            fragments.append(f'<span class="{css_class}">{escaped_base}</span>')
        else:
            fragments.append(escaped_base)
    return "".join(fragments)


def _render_reference_position_markers(
    aligned_reference: str,
    interval: int = 10,
) -> str:
    reference_position = 0
    markers: list[str] = []
    for alignment_index, base in enumerate(aligned_reference):
        if base == "-":
            continue
        reference_position += 1
        if reference_position % interval == 0:
            markers.append(
                f'<span class="position-marker" '
                f'style="left: calc(11ch + {alignment_index}ch)">'
                f"{reference_position}</span>"
            )
    return "".join(markers)


def _reference_tandem_classes(
    aligned_reference: str,
    itd: ITD,
) -> list[str | None]:
    classes: list[str | None] = []
    reference_position = -1
    for base in aligned_reference:
        css_class = None
        if base != "-":
            reference_position += 1
            if itd.copied_segment_start <= reference_position <= itd.copied_segment_end:
                css_class = "tandem-region"
        classes.append(css_class)
    return classes


def _alignment_difference_classes(
    alignment: Alignment,
    itd,
) -> list[str | None]:
    classes: list[str | None] = []
    ref_pos = -1
    insertion_offsets: dict[int, int] = {}
    expected_sequence = _expected_insertion_sequence(itd)
    prefix_length = len(itd.spacer_prefix)
    tandem_length = len(itd.copied_segment_sequence)

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        css_class = None
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-" and read_base != ref_base:
                css_class = "diff"
        elif read_base != "-":
            insertion_site = ref_pos
            offset = insertion_offsets.get(insertion_site, 0)
            if insertion_site == itd.insertion.start and offset < len(
                expected_sequence
            ):
                if offset < prefix_length or offset >= prefix_length + tandem_length:
                    css_class = "spacer-region"
                    if read_base != expected_sequence[offset]:
                        css_class = "spacer-region insert-mismatch"
                else:
                    css_class = "inserted-region"
                    expected_base = expected_sequence[offset]
                    if read_base != expected_base:
                        css_class = "inserted-region insert-mismatch"
            insertion_offsets[insertion_site] = offset + 1

        classes.append(css_class)

    return classes


def _alignment_comparison_classes(
    alignment,
    baseline_alignment,
) -> list[str | None]:
    baseline_reference_bases, baseline_insertions = _alignment_features(
        baseline_alignment.aligned_reference,
        baseline_alignment.aligned_read,
    )
    classes: list[str | None] = []
    ref_pos = -1
    insertion_offsets: dict[int, int] = {}

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        css_class = None
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-":
                baseline_base = baseline_reference_bases.get(ref_pos)
                if baseline_base is not None and baseline_base != read_base:
                    css_class = "diff"
        elif read_base != "-":
            insertion_site = ref_pos
            offset = insertion_offsets.get(insertion_site, 0)
            baseline_insertion = baseline_insertions.get(insertion_site, "")
            if offset < len(baseline_insertion):
                if baseline_insertion[offset] != read_base:
                    css_class = "diff"
            else:
                css_class = "insert"
            insertion_offsets[insertion_site] = offset + 1

        classes.append(css_class)

    return classes


def _alignment_features(
    aligned_reference: str,
    aligned_read: str,
) -> tuple[dict[int, str], dict[int, str]]:
    reference_bases: dict[int, str] = {}
    insertions: dict[int, list[str]] = {}
    ref_pos = -1

    for read_base, ref_base in zip(aligned_read, aligned_reference, strict=True):
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-":
                reference_bases[ref_pos] = read_base
            continue
        if read_base != "-":
            insertions.setdefault(ref_pos, []).append(read_base)

    return (
        reference_bases,
        {site: "".join(bases) for site, bases in insertions.items()},
    )


def _expected_insertion_sequence(itd: ITD) -> str:
    return f"{itd.spacer_prefix}{itd.copied_segment_sequence}{itd.spacer_suffix}"
