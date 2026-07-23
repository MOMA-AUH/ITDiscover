"""Command-line interface for ITDiscover."""

import argparse
import csv
import html
from pathlib import Path
from typing import TextIO

from . import __version__
from .alignment import align_read_to_reference
from .calls import (
    ITDCall,
    ITDFilter,
    UniqueSupportRepresentative,
    call_exact_itds_with_representatives,
    call_fuzzy_itds_with_representatives,
)
from .fastq import read_paired_fastq
from .insertions import Alignment
from .itds import ITD
from .reads import ReadTrimSettings, preprocess_fragments


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
        help="Optional forward primer sequence to trim from the 5' end of R1 reads.",
    )
    parser.add_argument(
        "--reverse-primer",
        help=(
            "Optional reverse primer sequence as it occurs at the 5' end of raw R2; "
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
        type=int,
        default=6,
        help="Minimum insertion length to consider.",
    )
    parser.add_argument(
        "--max-mismatches",
        type=_non_negative_int,
        help=(
            "Maximum mismatches allowed in the copied tandem tract; "
            "0 is equivalent to exact mode."
        ),
    )
    parser.add_argument(
        "--min-support-count",
        type=int,
        default=1,
        help="Minimum fragment support count required for a call to pass filtering.",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=0,
        help="Minimum coverage required for a call to pass filtering.",
    )
    parser.add_argument(
        "--min-vaf",
        type=float,
        default=0.0,
        help="Minimum VAF required for a call to pass filtering.",
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
    return _run_call_command(args)


def _run_call_command(args: argparse.Namespace) -> int:
    reference = _read_single_sequence_fasta(Path(args.reference))
    fragments = read_paired_fastq(args.r1, args.r2)
    trimming = _build_trim_settings(args)
    processed_reads = preprocess_fragments(
        fragments,
        min_length=args.min_read_length,
        min_mean_quality=args.min_mean_quality,
        trimming=trimming,
    )
    alignments = [
        align_read_to_reference(read, reference)
        for read in processed_reads
    ]
    filters = ITDFilter(
        min_support_count=args.min_support_count,
        min_coverage=args.min_coverage,
        min_vaf=args.min_vaf,
    )
    if args.max_mismatches is None:
        calls, representatives = call_exact_itds_with_representatives(
            alignments,
            reference,
            min_insert_length=args.min_insert_length,
            filters=filters,
        )
    else:
        calls, representatives = call_fuzzy_itds_with_representatives(
            alignments,
            reference,
            max_mismatches=args.max_mismatches,
            min_insert_length=args.min_insert_length,
            filters=filters,
        )
    if args.output:
        _write_unique_support_alignment_html_report(
            args.output,
            calls,
            representatives,
            filters=filters,
            max_mismatches=0 if args.max_mismatches is None else args.max_mismatches,
        )
    if args.output_tsv:
        _write_tsv_call_report(
            args.output_tsv,
            calls,
            max_mismatches=0 if args.max_mismatches is None else args.max_mismatches,
            min_support_count=filters.min_support_count,
            min_coverage=filters.min_coverage,
            min_vaf=filters.min_vaf,
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


def _build_trim_settings(args: argparse.Namespace) -> ReadTrimSettings | None:
    if not any(
        getattr(args, field) is not None
        for field in (
            "forward_primer",
            "reverse_primer",
        )
    ):
        return None
    return ReadTrimSettings(
        forward_primer=args.forward_primer,
        reverse_primer=args.reverse_primer,
    )


def _format_filter_reasons(call: ITDCall) -> str:
    return "." if not call.filter_reasons else ";".join(call.filter_reasons)


def _read_single_sequence_fasta(path: Path) -> str:
    with path.open(mode="rt", encoding="utf-8") as handle:
        sequences = list(_iter_fasta_sequences(handle))
    if not sequences:
        raise ValueError("reference FASTA does not contain a sequence")
    if len(sequences) > 1:
        raise ValueError("reference FASTA must contain exactly one sequence")
    return sequences[0]


def _iter_fasta_sequences(handle: TextIO) -> list[str]:
    sequences: list[str] = []
    current_parts: list[str] = []

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_parts:
                sequences.append("".join(current_parts))
                current_parts = []
            continue
        current_parts.append(line)

    if current_parts:
        sequences.append("".join(current_parts))
    return sequences


def _write_unique_support_alignment_html_report(
    path: Path,
    calls: list[ITDCall],
    representatives: list[UniqueSupportRepresentative],
    *,
    filters: ITDFilter | None = None,
    max_mismatches: int | None = None,
) -> None:
    representatives_by_key: dict[
        tuple[int, int, str, str, str, bool], list[UniqueSupportRepresentative]
    ] = {}
    for representative in representatives:
        key = (
            representative.itd.insertion.start,
            representative.itd.tandem_start,
            representative.itd.tandem_sequence,
            representative.itd.spacer_prefix,
            representative.itd.spacer_suffix,
            representative.itd.insertion.trailing,
        )
        representatives_by_key.setdefault(key, []).append(representative)

    sections: list[str] = []
    ordered_calls = sorted(
        calls,
        key=lambda call: (
            -call.support_count,
            call.itd.insertion.start,
            call.itd.tandem_start,
            call.itd.tandem_sequence,
            call.itd.spacer_prefix,
            call.itd.spacer_suffix,
            call.itd.insertion.sequence,
        ),
    )
    for call in ordered_calls:
        key = (
            call.itd.insertion.start,
            call.itd.tandem_start,
            call.itd.tandem_sequence,
            call.itd.spacer_prefix,
            call.itd.spacer_suffix,
            call.itd.insertion.trailing,
        )
        call_representatives = representatives_by_key.get(key, [])
        sections.append(_render_html_call_section(call, call_representatives))

    thresholds_section = _render_html_thresholds_section(
        filters=filters,
        max_mismatches=max_mismatches,
    )

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
  __THRESHOLDS__
  <div class="legend">
    <span class="legend-item"><span class="legend-chip tandem-region">T</span> tandem sequence</span>
    <span class="legend-item"><span class="legend-chip inserted-region">I</span> inserted sequence</span>
    <span class="legend-item"><span class="legend-chip spacer-region">S</span> spacer sequence</span>
    <span class="legend-item"><span class="legend-chip diff">A</span> mismatches</span>
  </div>
  __SECTIONS__
</body>
</html>
"""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        document.replace("__THRESHOLDS__", thresholds_section).replace(
            "__SECTIONS__",
            "\n".join(sections),
        ),
        encoding="utf-8",
    )


def _write_tsv_call_report(
    path: Path,
    calls: list[ITDCall],
    *,
    max_mismatches: int,
    min_support_count: int,
    min_coverage: int,
    min_vaf: float,
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
                "Insertion Start",
                "Tandem Start",
                "Tandem End",
                "Tandem Sequence",
                "Spacer Prefix",
                "Spacer Suffix",
                "Insertion Sequence",
                "Read-Edge Observation",
                "Support Count",
                "Coverage",
                "VAF",
                "Min Support Count",
                "Min Coverage",
                "Min VAF",
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
                    call.itd.tandem_start,
                    call.itd.tandem_end,
                    call.itd.tandem_sequence,
                    call.itd.spacer_prefix or "-",
                    call.itd.spacer_suffix or "-",
                    call.itd.insertion.sequence,
                    "Yes" if call.itd.is_partial_observation else "No",
                    call.support_count,
                    call.coverage,
                    f"{call.vaf:.6f}",
                    min_support_count,
                    min_coverage,
                    f"{min_vaf:.6f}",
                ]
            )


def _render_html_thresholds_section(
    *,
    filters: ITDFilter | None,
    max_mismatches: int,
) -> str:
    items: list[tuple[str, str]] = []
    if filters is not None:
        items.extend(
            [
                ("Min support count", str(filters.min_support_count)),
                ("Min coverage", str(filters.min_coverage)),
                ("Min VAF", f"{filters.min_vaf:.6f}"),
            ]
        )
    items.append(("Max mismatches", str(max_mismatches)))

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


def _render_html_call_section(
    call: ITDCall,
    representatives: list[UniqueSupportRepresentative],
) -> str:
    summary = (
        ('Insertion Start', str(call.itd.insertion.start)),
        ('Tandem Start', str(call.itd.tandem_start)),
        ('Sequence', call.itd.tandem_sequence),
        ('Spacer Prefix', call.itd.spacer_prefix or "-"),
        ('Spacer Suffix', call.itd.spacer_suffix or "-"),
        (
            'Read-Edge Observation',
            'Yes — partial; full ITD not reconstructed'
            if call.itd.is_partial_observation
            else 'No',
        ),
        ('Support Count', str(call.support_count)),
        ('Coverage', str(call.coverage)),
        ('VAF', f"{call.vaf:.6f}"),
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
    tandem_length = len(itd.tandem_sequence)

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
            if itd.tandem_start <= ref_pos <= itd.tandem_end:
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
    tandem_length = len(itd.tandem_sequence)

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
    return f"{itd.spacer_prefix}{itd.tandem_sequence}{itd.spacer_suffix}"
