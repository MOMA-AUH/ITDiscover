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
    "Reference-local, zero-based. Insertion coordinate is the reference base "
    "immediately before the insertion (-1 means before the first base). Copied "
    "segment start is zero-based; copied segment end is zero-based and inclusive. "
    "The copied segment is immediately before the insertion when it ends at the "
    "insertion coordinate, and immediately after when it starts at the following "
    "reference base. Before/after describes this coordinate representation, not "
    "a biological direction of copying."
)


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
        "--max-copy-mismatches",
        dest="max_copy_mismatches",
        type=_non_negative_int,
        help=(
            "Maximum mismatches allowed in the copied reference segment; "
            "0 is equivalent to exact mode. ITD detection happens before any "
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
        "--consolidation-max-allele-mismatches",
        dest="consolidation_max_allele_mismatches",
        type=_non_negative_int,
        default=1,
        help=(
            "Maximum positional mismatches between the complete ALT sequences "
            "of an already-detected minor allele and its anchor. This advanced "
            "setting does not affect whether an insertion qualifies as an ITD."
        ),
    )
    consolidation_group.add_argument(
        "--consolidation-max-breakpoint-shift",
        type=_non_negative_int,
        default=6,
        help="Maximum reference-base shift between minor and anchor breakpoints.",
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
        "--min-junction-quality",
        type=_non_negative_int,
        default=30,
        help="Minimum Phred quality for inserted bases and junction anchors.",
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
        help="Optional path for an HTML report with one representative alignment per unique support pattern.",
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
        min_junction_quality=args.min_junction_quality,
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
    if args.max_copy_mismatches is None:
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
            max_mismatches=args.max_copy_mismatches,
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
        _write_unique_support_alignment_html_report(
            args.output,
            calls,
            representatives,
            filters=filters,
            max_mismatches=(
                0
                if args.max_copy_mismatches is None
                else args.max_copy_mismatches
            ),
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
    if args.output_tsv:
        _write_tsv_call_report(
            args.output_tsv,
            calls,
            max_mismatches=(
                0
                if args.max_copy_mismatches is None
                else args.max_copy_mismatches
            ),
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
        max_allele_mismatches=args.consolidation_max_allele_mismatches,
        max_breakpoint_shift=args.consolidation_max_breakpoint_shift,
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
            _write_unique_support_alignment_html_report(
                args.output,
                [],
                [],
                max_mismatches=(
                    0
                    if args.max_copy_mismatches is None
                    else args.max_copy_mismatches
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
    if args.output_tsv:
        try:
            _write_tsv_call_report(
                args.output_tsv,
                [],
                max_mismatches=(
                    0
                    if args.max_copy_mismatches is None
                    else args.max_copy_mismatches
                ),
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
            f"breakpoint_shift={member.breakpoint_shift},"
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


def _write_unique_support_alignment_html_report(
    path: Path,
    calls: list[ITDCall],
    representatives: list[UniqueSupportRepresentative],
    *,
    filters: ITDFilter | None = None,
    max_mismatches: int | None = None,
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
    representatives_by_key: dict[object, list[UniqueSupportRepresentative]] = {}
    for representative in representatives:
        representatives_by_key.setdefault(
            representative.canonical_allele,
            [],
        ).append(representative)

    sections: list[str] = []
    ordered_calls = sorted(
        calls,
        key=lambda call: (
            -call.mutant_fragment_count,
            call.itd.insertion.start,
            call.itd.copied_segment_start,
            call.itd.copied_segment_sequence,
            call.itd.spacer_prefix,
            call.itd.spacer_suffix,
            call.itd.insertion.sequence,
        ),
    )
    for call in ordered_calls:
        call_representatives = representatives_by_key.get(
            call.canonical_allele,
            [],
        )
        sections.append(_render_html_call_section(call, call_representatives))

    thresholds_section = _render_html_thresholds_section(
        filters=filters,
        max_mismatches=max_mismatches,
        alignment_filters=alignment_filters,
        insertion_filters=insertion_filters,
        min_insert_length=min_insert_length,
        min_copied_segment_length=min_copied_segment_length,
        require_in_frame=require_in_frame,
        consolidation=consolidation,
    )
    sample_summary = _render_html_sample_summary(sample_result, qc_thresholds)
    reference_summary = _render_html_reference_summary(
        reference_id,
        reference_length,
        reference_sha256,
    )
    empty_state = "" if sections else _render_html_empty_state(sample_result)

    document = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>ITDiscover Report</title>
  <style>
    :root {
      color-scheme: light;
      --ink: #18232f;
      --muted: #5a6875;
      --line: #cfd8df;
      --panel: #f7fafc;
      --tandem-bg: #dbeafe;
      --tandem-fg: #12315d;
      --inserted-bg: #dbeafe;
      --inserted-fg: #12315d;
      --spacer-bg: #fef3c7;
      --spacer-fg: #92400e;
      --mismatch-bg: #fee2e2;
      --mismatch-fg: #b91c1c;
    }
    body {
      margin: 24px;
      color: var(--ink);
      font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
      line-height: 1.4;
    }
    h1, h2, h3 {
      margin: 0 0 12px;
    }
    .itd {
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 16px;
      margin-bottom: 20px;
      background: white;
    }
    .summary {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
      gap: 10px 16px;
      margin-bottom: 16px;
    }
    .summary div {
      background: var(--panel);
      border-radius: 4px;
      padding: 8px 10px;
    }
    .summary dt {
      margin: 0 0 4px;
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }
    .summary dd {
      margin: 0;
      font-size: 16px;
      font-weight: 600;
    }
    .support {
      border-top: 1px solid var(--line);
      padding-top: 14px;
      margin-top: 14px;
    }
    .support-header {
      display: flex;
      gap: 16px;
      flex-wrap: wrap;
      margin-bottom: 10px;
      align-items: baseline;
    }
    .support-header strong {
      font-size: 16px;
    }
    .support-meta {
      color: var(--muted);
      font-size: 14px;
    }
    .legend {
      display: flex;
      gap: 12px 18px;
      flex-wrap: wrap;
      margin: 0 0 20px;
      color: var(--ink);
      font-size: 15px;
    }
    .thresholds {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--panel);
      padding: 12px 14px;
      margin: 0 0 20px;
    }
    .thresholds-title {
      margin: 0 0 10px;
      color: var(--ink);
      font-size: 14px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.06em;
    }
    .thresholds-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 8px 14px;
    }
    .thresholds-item {
      display: grid;
      gap: 2px;
    }
    .thresholds-label {
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }
    .thresholds-value {
      font-size: 15px;
      font-weight: 600;
    }
    .quantification-note {
      color: var(--muted);
      font-size: 14px;
      line-height: 1.4;
      margin: -4px 0 20px;
    }
    .sample-result, .empty-state {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--panel);
      padding: 12px 14px;
      margin: 0 0 20px;
    }
    .sample-result h2 { margin-top: 0; }
    .legend-item {
      display: inline-flex;
      gap: 8px;
      align-items: center;
    }
    .legend-chip {
      min-width: 1.6em;
      border: 1px solid rgb(24 35 47 / 18%);
      border-radius: 4px;
      padding: 2px 7px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 16px;
      line-height: 1.25;
      text-align: center;
    }
    .representative-title,
    .pileup-title {
      margin: 14px 0 8px;
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }
    .signature {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: var(--panel);
      border-radius: 4px;
      padding: 6px 8px;
      margin-bottom: 10px;
      display: inline-block;
    }
    .alignment-block {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      overflow-x: auto;
      background: #fbfdff;
      border: 1px solid var(--line);
      border-radius: 4px;
      padding: 10px;
    }
    .alignment-row {
      white-space: pre;
    }
    .match-comparison {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      background: var(--panel);
      border-radius: 4px;
      padding: 8px 10px;
      margin: 0 0 10px;
      overflow-x: auto;
    }
    .match-row {
      white-space: pre;
    }
    .pileup {
      margin-top: 12px;
    }
    .pileup-table {
      border-collapse: collapse;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      min-width: min(100%, 560px);
    }
    .pileup-table th,
    .pileup-table td {
      border-bottom: 1px solid var(--line);
      padding: 5px 8px;
      text-align: left;
      white-space: pre;
    }
    .pileup-table th {
      color: var(--muted);
      font-weight: 600;
    }
    .label {
      color: var(--muted);
    }
    .diff {
      background: var(--mismatch-bg);
      color: var(--mismatch-fg);
      font-weight: 700;
    }
    .tandem-region {
      background: var(--tandem-bg);
      color: var(--tandem-fg);
      font-weight: 700;
    }
    .inserted-region {
      background: var(--inserted-bg);
      color: var(--inserted-fg);
      font-weight: 700;
    }
    .spacer-region {
      background: var(--spacer-bg);
      color: var(--spacer-fg);
      font-weight: 700;
    }
    .insert-mismatch {
      background: var(--mismatch-bg);
      color: var(--mismatch-fg);
      font-weight: 700;
    }
  </style>
</head>
<body>
  <h1>ITDiscover Report</h1>
  __REFERENCE_SUMMARY__
  __SAMPLE_SUMMARY__
  __THRESHOLDS__
  <p class="quantification-note">Observed mutant-fragment fraction = mutant
  fragments / (mutant + wild-type fragments). Both allele states require the
  configured high-quality junction anchors. Conflicting, unresolved, and
  not-informative fragments are reported separately. Overlapping mates count
  once per fragment. PCR duplicates are not collapsed unless they already
  share a fragment ID. Direction bias compares the mutant fraction among
  callable junction opportunities in R1 with the corresponding fraction in
  R2; it is not evaluated unless both directions meet the configured
  opportunity threshold.</p>
  <div class="legend">
    <span class="legend-item"><span class="legend-chip tandem-region">C</span> copied reference segment</span>
    <span class="legend-item"><span class="legend-chip inserted-region">I</span> inserted sequence</span>
    <span class="legend-item"><span class="legend-chip spacer-region">S</span> spacer sequence</span>
    <span class="legend-item"><span class="legend-chip diff">A</span> mismatches</span>
  </div>
  __EMPTY_STATE__
  __SECTIONS__
</body>
</html>
"""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        document.replace("__REFERENCE_SUMMARY__", reference_summary)
        .replace("__SAMPLE_SUMMARY__", sample_summary)
        .replace("__THRESHOLDS__", thresholds_section)
        .replace("__EMPTY_STATE__", empty_state)
        .replace("__SECTIONS__", "\n".join(sections)),
        encoding="utf-8",
    )


def _write_tsv_call_report(
    path: Path,
    calls: list[ITDCall],
    *,
    max_mismatches: int,
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
                "Max Mismatches",
                "Insertion After Reference Base (0-based; -1=before first)",
                "Copied Segment Start (0-based)",
                "Copied Segment End (0-based, inclusive)",
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
                "Min Junction Quality",
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
                "Consolidation Max Allele Mismatches",
                "Consolidation Max Breakpoint Shift",
                "Consolidation Max Minor/Anchor Support Ratio",
                "Consolidation Min Anchor Fragment Count",
                "Consolidated Minor Allele Count",
                "Consolidated Minor Raw Fragment Support",
                "Consolidated Minor Alleles",
            ]
        )
        mode = "exact" if max_mismatches == 0 else "fuzzy"
        for index, call in enumerate(calls, start=1):
            writer.writerow(
                [
                    index,
                    call.status,
                    _format_filter_reasons(call),
                    mode,
                    max_mismatches,
                    call.itd.insertion.start,
                    call.itd.copied_segment_start,
                    call.itd.copied_segment_end,
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
                        insertion_filters.min_junction_quality
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
                    consolidation.max_allele_mismatches,
                    consolidation.max_breakpoint_shift,
                    f"{consolidation.max_minor_to_anchor_support_ratio:.6f}",
                    consolidation.min_anchor_fragment_count,
                    len(call.consolidated_members),
                    call.consolidated_minor_fragment_count,
                    _format_consolidated_members(call),
                ]
            )
        if not calls:
            empty_call_values: list[object] = ["."] * 34
            empty_call_values[3] = mode
            empty_call_values[4] = max_mismatches
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
                empty_call_values[32] = insertion_filters.min_junction_quality
                empty_call_values[33] = insertion_filters.junction_flank_size
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
                    consolidation.max_allele_mismatches,
                    consolidation.max_breakpoint_shift,
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
            "No ITD candidates were called, but the result is indeterminate "
            "because sample QC did not pass."
        )
    else:
        message = "No ITD candidates were called."
    return f'<section class="empty-state">{html.escape(message)}</section>'


def _render_html_sample_summary(
    result: SampleResult | None,
    thresholds: SampleQCThresholds | None,
) -> str:
    if result is None:
        return ""
    values: list[tuple[str, str]] = [
        ("Sample", result.sample_id),
        ("Analysis Status", result.analysis_status),
        ("QC Status", result.qc_status),
        ("Outcome", result.outcome),
        ("QC Reasons", "; ".join(result.qc_reasons) or "None"),
        ("Passing Calls", str(result.passing_call_count)),
        ("Filtered Candidates", str(result.filtered_candidate_count)),
    ]
    if result.error_message:
        values.append(("Analysis Error", result.error_message))
    if result.preprocessing is not None:
        metrics = result.preprocessing
        values.extend(
            [
                ("Input Fragments", str(metrics.input_fragment_count)),
                ("Input Reads", str(metrics.input_read_count)),
                ("Usable Fragments", str(metrics.usable_fragment_count)),
                (
                    "Preprocessing-passing Reads (R1/R2)",
                    f"{metrics.passing_read_count} "
                    f"({metrics.passing_forward_read_count}/"
                    f"{metrics.passing_reverse_read_count})",
                ),
                ("Primer-failed Reads", str(metrics.primer_failed_read_count)),
                ("Length-failed Reads", str(metrics.length_failed_read_count)),
                ("Quality-failed Reads", str(metrics.quality_failed_read_count)),
            ]
        )
        if metrics.primer_retained_forward_reads is not None:
            values.append(
                (
                    "Forward Primer-retained Reads",
                    str(metrics.primer_retained_forward_reads),
                )
            )
        if metrics.primer_retained_reverse_reads is not None:
            values.append(
                (
                    "Reverse Primer-retained Reads",
                    str(metrics.primer_retained_reverse_reads),
                )
            )
    if result.alignment is not None:
        metrics = result.alignment
        values.extend(
            [
                ("Alignment-passing Reads", str(metrics.passing_read_count)),
                ("Alignment Pass Fraction", f"{metrics.pass_fraction:.1%}"),
                (
                    "Alignment-passing Reads (R1/R2)",
                    f"{metrics.passing_forward_read_count}/"
                    f"{metrics.passing_reverse_read_count}",
                ),
            ]
        )
    if result.coverage is not None:
        values.append(
            (
                "Inter-base Coverage (min/median/max)",
                f"{result.coverage.minimum}/"
                f"{result.coverage.median:.1f}/"
                f"{result.coverage.maximum}",
            )
        )
    if thresholds is not None:
        values.extend(
            [
                ("QC Min Usable Fragments", str(thresholds.min_usable_fragment_count)),
                (
                    "QC Min Reads per Direction",
                    str(thresholds.min_passing_reads_per_direction),
                ),
                (
                    "QC Min Alignment Pass Fraction",
                    f"{thresholds.min_alignment_pass_fraction:.1%}",
                ),
                (
                    "QC Min Median Coverage",
                    str(thresholds.min_median_interbase_coverage),
                ),
                (
                    "QC Min Primer Retention",
                    f"{thresholds.min_primer_retention_fraction:.1%}",
                ),
            ]
        )
    value_html = "".join(
        f"<div><dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd></div>"
        for label, value in values
    )
    return (
        '<section class="sample-result">'
        '<h2>Sample Result and QC</h2>'
        f'<dl class="summary">{value_html}</dl>'
        "</section>"
    )


def _render_html_thresholds_section(
    *,
    filters: ITDFilter | None,
    max_mismatches: int,
    alignment_filters: AlignmentEvidenceFilter | None = None,
    insertion_filters: InsertionEvidenceFilter | None = None,
    min_insert_length: int = 6,
    min_copied_segment_length: int = 6,
    require_in_frame: bool = True,
    consolidation: ITDConsolidationSettings = ITDConsolidationSettings(),
) -> str:
    items: list[tuple[str, str]] = []
    if filters is not None:
        items.extend(
            [
                (
                    "Min mutant fragments",
                    str(filters.min_mutant_fragment_count),
                ),
                (
                    "Min informative fragments",
                    str(filters.min_informative_fragment_count),
                ),
                (
                    "Min mutant-fragment fraction",
                    f"{filters.min_observed_mutant_fragment_fraction:.3%}",
                ),
                (
                    "Max directional mutant-fraction share",
                    f"{filters.max_directional_mutant_fraction_share:.3f}",
                ),
                (
                    "Min opportunities per direction",
                    str(filters.min_directional_opportunities),
                ),
            ]
        )
    items.extend(
        [
            ("Max copied-segment mismatches", str(max_mismatches)),
            ("Min insert length", str(min_insert_length)),
            ("Min copied-segment length", str(min_copied_segment_length)),
            ("Require in-frame insertions", "Yes" if require_in_frame else "No"),
            (
                "Minor-variant consolidation",
                "Enabled" if consolidation.enabled else "Disabled",
            ),
            (
                "Consolidation max allele mismatches",
                str(consolidation.max_allele_mismatches),
            ),
            (
                "Consolidation max breakpoint shift",
                str(consolidation.max_breakpoint_shift),
            ),
            (
                "Consolidation max minor/anchor support ratio",
                f"{consolidation.max_minor_to_anchor_support_ratio:.3f}",
            ),
            (
                "Consolidation min anchor fragments",
                str(consolidation.min_anchor_fragment_count),
            ),
        ]
    )
    if alignment_filters is not None:
        items.extend(
            [
                ("Min alignment identity", f"{alignment_filters.min_identity:.3f}"),
                (
                    "Min on-target fraction",
                    f"{alignment_filters.min_on_target_fraction:.3f}",
                ),
                (
                    "Min alignment score",
                    str(alignment_filters.min_score)
                    if alignment_filters.min_score is not None
                    else "Not set",
                ),
                (
                    "Reject ambiguous alignments",
                    "Yes" if alignment_filters.reject_ambiguous else "No",
                ),
            ]
        )
    if insertion_filters is not None:
        items.extend(
            [
                (
                    "Min junction quality",
                    str(insertion_filters.min_junction_quality),
                ),
                ("Junction flank size", str(insertion_filters.junction_flank_size)),
            ]
        )

    item_html = "".join(
        (
            '<div class="thresholds-item">'
            f'<span class="thresholds-label">{html.escape(label)}</span>'
            f'<span class="thresholds-value">{html.escape(value)}</span>'
            "</div>"
        )
        for label, value in items
    )
    return (
        '<section class="thresholds">'
        '<div class="thresholds-title"><strong>CALL THRESHOLDS</strong></div>'
        f'<div class="thresholds-grid">{item_html}</div>'
        '</section>'
    )


def _render_html_reference_summary(
    reference_id: str | None,
    reference_length: int | None,
    reference_sha256: str | None,
) -> str:
    values = (
        ("Reference FASTA Header", reference_id or "Not provided"),
        (
            "Reference Length",
            str(reference_length) if reference_length is not None else "Not provided",
        ),
        ("Reference Sequence SHA-256", reference_sha256 or "Not provided"),
        ("Coordinate Convention", COORDINATE_CONVENTION),
    )
    value_html = "".join(
        f"<div><dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd></div>"
        for label, value in values
    )
    return (
        '<section class="sample-result">'
        '<h2>Reference and Coordinates</h2>'
        f'<dl class="summary">{value_html}</dl>'
        "</section>"
    )


def _render_html_call_section(
    call: ITDCall,
    representatives: list[UniqueSupportRepresentative],
) -> str:
    summary = (
        (
            'Insertion After Reference Base (0-based)',
            str(call.itd.insertion.start),
        ),
        ('Copied Segment Start (0-based)', str(call.itd.copied_segment_start)),
        (
            'Copied Segment End (0-based, inclusive)',
            str(call.itd.copied_segment_end),
        ),
        (
            'Copied Segment Location',
            f"Immediately {call.itd.copied_segment_location} the insertion",
        ),
        ('Copied Segment Sequence', call.itd.copied_segment_sequence),
        ('Spacer Prefix', call.itd.spacer_prefix or "-"),
        ('Spacer Suffix', call.itd.spacer_suffix or "-"),
        (
            'Read-Edge Observation',
            'Yes — partial; full ITD not reconstructed'
            if call.itd.is_partial_observation
            else 'No',
        ),
        ('Mutant Fragments', str(call.mutant_fragment_count)),
        ('Consolidated Minor Alleles', str(len(call.consolidated_members))),
        (
            'Consolidated Minor Raw Fragment Support',
            str(call.consolidated_minor_fragment_count),
        ),
        ('Consolidation Audit', _format_consolidated_members(call)),
        ('Wild-type Fragments', str(call.wild_type_fragment_count)),
        ('Informative Fragments', str(call.informative_fragment_count)),
        (
            'R1 Evidence',
            _format_directional_evidence(
                call.r1_mutant_count,
                call.r1_opportunity_count,
            ),
        ),
        (
            'R2 Evidence',
            _format_directional_evidence(
                call.r2_mutant_count,
                call.r2_opportunity_count,
            ),
        ),
        ('Concordant Fragments', str(call.concordant_fragment_count)),
        ('Single-mate Fragments', str(call.single_mate_fragment_count)),
        (
            'Conflicting Fragments',
            str(call.conflicting_fragment_count),
        ),
        (
            'Unresolved Fragments',
            str(call.unresolved_fragment_count),
        ),
        (
            'Not-informative Fragments',
            str(call.not_informative_fragment_count),
        ),
        (
            'Observed Mutant-fragment Fraction',
            (
                f"{call.observed_mutant_fragment_fraction:.1%} "
                f"({call.mutant_fragment_count}/"
                f"{call.informative_fragment_count} informative fragments)"
            ),
        ),
        ('Status', call.status),
        ('Filter Reasons', _format_filter_reasons(call)),
    )
    summary_html = "".join(
        f"<div><dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd></div>"
        for label, value in summary
    )
    representative = _best_representative(representatives)
    support_html = (
        _render_html_support_block(representative)
        if representative is not None
        else ""
    )
    return (
        f'<section class="itd">'
        f'<dl class="summary">{summary_html}</dl>'
        f"{support_html}"
        f"</section>"
    )


def _best_representative(
    representatives: list[UniqueSupportRepresentative],
) -> UniqueSupportRepresentative | None:
    if not representatives:
        return None
    return min(
        representatives,
        key=lambda representative: (
            representative.mismatches,
            -representative.support_count,
            representative.signature,
            representative.alignment.read_id,
        ),
    )


def _render_html_support_block(
    representative: UniqueSupportRepresentative,
) -> str:
    alignment = representative.alignment
    reference_classes = _reference_tandem_classes(
        alignment.aligned_reference,
        representative.itd,
    )
    reference_html = _highlight_alignment_differences(
        alignment.aligned_reference,
        comparison_classes=reference_classes,
    )
    comparison_classes = _alignment_difference_classes(
        alignment,
        representative.itd,
    )
    read_html = _highlight_alignment_differences(
        alignment.aligned_read,
        comparison_classes=comparison_classes,
    )
    return (
        '<section class="support">'
        '<div class="representative-title">Representative alignment</div>'
        f'<div class="support-header"><strong>{html.escape(alignment.read_id)}</strong></div>'
        '<div class="alignment-block">'
        f'<div class="alignment-row"><span class="label">reference  </span>{reference_html}</div>'
        f'<div class="alignment-row"><span class="label">read       </span>{read_html}</div>'
        '</div>'
        f'{_render_insert_sequence_pileup(representative)}'
        '</section>'
    )


def _render_insert_sequence_pileup(
    representative: UniqueSupportRepresentative,
) -> str:
    if not representative.insert_sequence_supports:
        return ""

    rows = "".join(
        (
            "<tr>"
            f"<td>{_highlight_inserted_sequence(support.sequence, representative.itd)}</td>"
            f"<td>{support.mismatches}</td>"
            f"<td>{support.support_count}</td>"
            "</tr>"
        )
        for support in representative.insert_sequence_supports
    )
    return (
        '<div class="pileup">'
        '<div class="pileup-title">Inserted sequence pileup</div>'
        '<table class="pileup-table">'
        '<thead><tr><th>Inserted sequence</th><th>Mismatches</th><th>Count</th></tr></thead>'
        f"<tbody>{rows}</tbody>"
        "</table>"
        '</div>'
    )


def _highlight_sequence_mismatches(observed: str, expected: str) -> str:
    fragments: list[str] = []
    for observed_base, expected_base in zip(observed, expected, strict=True):
        escaped_base = html.escape(observed_base)
        if observed_base == expected_base:
            fragments.append(escaped_base)
            continue
        fragments.append(f'<span class="insert-mismatch">{escaped_base}</span>')
    return "".join(fragments)


def _highlight_inserted_sequence(sequence: str, itd: ITD) -> str:
    fragments: list[str] = []
    expected_sequence = _expected_insertion_sequence(itd)
    prefix_length = len(itd.spacer_prefix)
    tandem_length = len(itd.copied_segment_sequence)

    for index, base in enumerate(sequence):
        escaped_base = html.escape(base)
        if index < prefix_length or index >= prefix_length + tandem_length:
            css_class = "spacer-region"
            if base != expected_sequence[index]:
                css_class = "spacer-region insert-mismatch"
            fragments.append(f'<span class="{css_class}">{escaped_base}</span>')
            continue

        css_class = "inserted-region"
        if base != expected_sequence[index]:
            css_class = "inserted-region insert-mismatch"
        fragments.append(f'<span class="{css_class}">{escaped_base}</span>')

    return "".join(fragments)


def _highlight_alignment_differences(
    sequence: str,
    *,
    comparison_classes: list[str | None] | None,
) -> str:
    fragments: list[str] = []
    for index, base in enumerate(sequence):
        escaped_base = html.escape(base)
        css_class = None
        if comparison_classes is not None and index < len(comparison_classes):
            css_class = comparison_classes[index]
        if css_class is not None:
            fragments.append(f'<span class="{css_class}">{escaped_base}</span>')
            continue
        fragments.append(escaped_base)
    return "".join(fragments)


def _reference_tandem_classes(
    aligned_reference: str,
    itd,
) -> list[str | None]:
    classes: list[str | None] = []
    ref_pos = -1

    for ref_base in aligned_reference:
        css_class = None
        if ref_base != "-":
            ref_pos += 1
            if itd.copied_segment_start <= ref_pos <= itd.copied_segment_end:
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
