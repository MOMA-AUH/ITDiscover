import csv
import hashlib
from pathlib import Path

import pytest

import itdiscover
import itdiscover.cli as cli
from itdiscover.alleles import CanonicalInsertionAllele
from itdiscover.calls import (
    ConsolidatedAlleleMember,
    ITDCall,
    ITDConsolidationSettings,
    UniqueSupportRepresentative,
)
from itdiscover.insertions import Alignment
from itdiscover.insertions import Insertion
from itdiscover.itds import ITD
from itdiscover.sequences import reverse_complement


def test_main_version(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == f"itdiscover {itdiscover.__version__}\n"


def test_main_requires_arguments(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main([])

    assert exc_info.value.code == 2
    assert "required" in capsys.readouterr().err


@pytest.mark.parametrize("missing_primer", ["--forward-primer", "--reverse-primer"])
def test_call_command_requires_both_primers(missing_primer, capsys) -> None:
    arguments = [
        "--reference",
        "reference.fasta",
        "--r1",
        "r1.fastq",
        "--r2",
        "r2.fastq",
        "--forward-primer",
        "AAA",
        "--reverse-primer",
        "TTT",
    ]
    primer_index = arguments.index(missing_primer)
    del arguments[primer_index : primer_index + 2]

    with pytest.raises(SystemExit) as exc_info:
        cli.main(arguments)

    assert exc_info.value.code == 2
    assert missing_primer in capsys.readouterr().err


def test_reversed_primer_exits_cleanly_with_correction(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(">reference\nAAACCCGGGTTT\n", encoding="utf-8")
    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        "@fragment/1\nGGGAAACCC\n+\nIIIIIIIII\n",
        encoding="utf-8",
    )
    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        "@fragment/2\nTTTAAACCCTTT\n+\nIIIIIIIIIIII\n",
        encoding="utf-8",
    )
    report_path = tmp_path / "report.tsv"

    assert cli.main(
        [
            "--reference",
            str(reference_path),
            "--r1",
            str(r1_path),
            "--r2",
            str(r2_path),
            "--forward-primer",
            "GGG",
            "--reverse-primer",
            "CCCAAA",
            "--output-tsv",
            str(report_path),
        ]
    ) == 2

    error = capsys.readouterr().err
    assert "Traceback" not in error
    assert "--reverse-primer matches no raw R2 reads" in error
    assert "'AAACCC' matches 1" in error
    rows = list(csv.DictReader(report_path.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    assert rows[0]["Analysis Error"].startswith("--reverse-primer matches no raw R2 reads")


def test_event_length_thresholds_must_be_positive(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "--reference",
                "reference.fasta",
                "--r1",
                "r1.fastq",
                "--r2",
                "r2.fastq",
                "--min-copied-segment-length",
                "0",
            ]
        )

    assert exc_info.value.code == 2
    assert "at least 1" in capsys.readouterr().err


def test_direction_and_copied_segment_options_use_clear_destinations() -> None:
    parser = cli.build_parser()
    args = parser.parse_args(
        [
            "--reference",
            "reference.fasta",
            "--r1",
            "r1.fastq",
            "--r2",
            "r2.fastq",
            "--forward-primer",
            "AAA",
            "--reverse-primer",
            "AAA",
            "--max-directional-mutant-fraction-share",
            "0.8",
            "--min-directional-opportunities",
            "7",
            "--min-copied-segment-length",
            "9",
            "--max-copy-mismatch-rate",
            "0.125",
            "--min-junction-anchor-quality",
            "31",
            "--min-insert-mean-quality",
            "29.5",
            "--min-insert-base-quality",
            "12",
            "--consolidate-minor-itd-variants",
            "--consolidation-max-allele-mismatch-rate",
            "0.2",
            "--consolidation-max-breakpoint-shift-rate",
            "0.5",
            "--consolidation-max-minor-support-ratio",
            "0.02",
            "--consolidation-min-anchor-fragment-count",
            "8",
        ]
    )

    assert args.max_directional_mutant_fraction_share == 0.8
    assert args.min_directional_opportunities == 7
    assert args.min_copied_segment_length == 9
    assert args.max_copy_mismatch_rate == 0.125
    assert args.min_junction_anchor_quality == 31
    assert args.min_insert_mean_quality == 29.5
    assert args.min_insert_base_quality == 12
    assert args.consolidate_minor_itd_variants is True
    assert args.consolidation_max_allele_mismatch_rate == 0.2
    assert args.consolidation_max_breakpoint_shift_rate == 0.5
    assert args.consolidation_max_minor_support_ratio == 0.02
    assert args.consolidation_min_anchor_fragment_count == 8
    help_text = parser.format_help()
    assert "advanced minor-allele consolidation:" in help_text
    assert "--max-copy-mismatch-rate" in help_text
    assert "--consolidation-max-allele-mismatch-rate" in help_text
    assert "--consolidation-max-breakpoint-shift-rate" in help_text


@pytest.mark.parametrize(
    "removed_option",
    [
        "--max-copy-mismatches",
        "--consolidation-max-allele-mismatches",
        "--consolidation-max-breakpoint-shift",
        "--min-junction-quality",
    ],
)
def test_removed_option_aliases_are_rejected(
    removed_option,
    capsys,
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.build_parser().parse_args(
            [
                "--reference",
                "reference.fasta",
                "--r1",
                "r1.fastq",
                "--r2",
                "r2.fastq",
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                removed_option,
                "1",
            ]
        )

    assert exc_info.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err


def test_documented_flt3_example_has_expected_interpretation(
    tmp_path,
    capsys,
) -> None:
    data_dir = Path(__file__).parent / "data" / "synthetic_flt3"
    html_path = tmp_path / "report.html"
    tsv_path = tmp_path / "calls.tsv"

    assert cli.main(
        [
            "--reference",
            str(data_dir / "reference.fa"),
            "--r1",
            str(data_dir / "synthetic_R1.fastq"),
            "--r2",
            str(data_dir / "synthetic_R2.fastq"),
            "--forward-primer",
            "GCAATTTAGGTATGAAAGCCAGC",
            "--reverse-primer",
            "CTTTCAGCATTTTGACGGCAACC",
            "--sample-id",
            "synthetic-flt3",
            "--output",
            str(html_path),
            "--output-tsv",
            str(tsv_path),
        ]
    ) == 0

    assert capsys.readouterr().out == ""
    rows = list(
        csv.DictReader(
            tsv_path.read_text(encoding="utf-8").splitlines(),
            delimiter="\t",
        )
    )
    assert len(rows) == 2
    passing, filtered = rows
    assert passing["Status"] == "PASS"
    assert passing["Sample ID"] == "synthetic-flt3"
    assert passing["Reference Length"] == "329"
    assert passing["Outcome"] == "ITD detected"
    assert passing["QC Status"] == "pass"
    assert passing["Insertion Sequence"] == "AGAGAATATGAATAT"
    insertion_coordinate = "Insertion After Reference Base"
    assert passing[insertion_coordinate] == "79"
    assert passing["Copied Segment Start"] == "80"
    assert passing["Copied Segment End"] == "94"
    assert passing["Copied Segment Location"] == "immediately after insertion"
    assert passing["Mutant Fragment Count"] == "3"
    assert passing["Wild-type Fragment Count"] == "9"
    assert passing["Informative Fragment Count"] == "12"
    assert passing["Observed Mutant-fragment Fraction"] == "0.250000"
    assert filtered["Status"] == "FAIL"
    assert filtered["Filter Reasons"] == "LOW_MUTANT_FRAGMENT_COUNT"

    html_report = html_path.read_text(encoding="utf-8")
    assert '<div class="sample-name">ITDiscover Report</div>' in html_report
    assert '<h2 class="outcome">synthetic-flt3</h2>' in html_report
    assert "ITD detected" not in html_report
    assert "QC Status:</span><span class=\"status-value status-value--success\">PASS" in html_report
    assert "AGAGAATATGAATAT" in html_report
    assert "Reference coordinates are 1-based" in html_report
    assert "80–94" in html_report
    assert "Copied segment location" not in html_report
    assert "Observed mutant-fragment fraction" in html_report
    assert "25.0%" in html_report
    assert "Mutant fragments</dt><dd>3" in html_report
    assert "Informative fragments</dt><dd>12" in html_report
    assert "Representative alignment" not in html_report
    assert '<section class="alignment">' in html_report
    assert '<details class="filtered-variants">' in html_report
    assert "Filtered variant (1)" in html_report
    assert "LOW-SUPPORT" in html_report
    assert 'class="filter-reason" title="Mutant-supporting fragment count is below the configured minimum."' in html_report
    assert "text-transform: uppercase" not in html_report
    assert 'class="alignment-ruler"' in html_report
    assert "copied reference</span>" in html_report
    assert "spacer</span>" in html_report
    assert "mismatch</span>" in html_report
    assert 'class="inserted-region"' in html_report
    assert "1 filtered candidate" in html_report
    assert "Inserted sequence pileup" not in html_report
    assert "CALL THRESHOLDS" not in html_report


def test_call_command_reports_exact_itd_from_paired_fastq(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            f"{reverse_complement('AAACCCGGGCCCGGGTTT')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""


def test_call_command_reports_fuzzy_itd_from_paired_fastq(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")
    report_path = tmp_path / "fuzzy-report.html"

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAAAAACCCGGGCCCGGATTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "AAAAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            f"AAA{reverse_complement('AAACCCGGGCCCGGATTT')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "AAAAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
                "--max-copy-mismatch-rate",
                str(1 / 6),
                "--output",
                str(report_path),
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""
    report = report_path.read_text(encoding="utf-8")
    assert "<title>ITDiscover Report</title>" in report
    assert "<h1>ITDiscover Report</h1>" not in report
    assert '<h2 class="outcome">' not in report
    assert "QC Status:</span><span class=\"status-value status-value--failure\">FAIL" in report
    assert "No passing ITD result is available" in report
    assert "See the TSV output for full QC metrics, thresholds, and audit details" in report
    assert "CALL THRESHOLDS" not in report
    assert "Representative alignment" not in report
    assert "Inserted sequence pileup" not in report
    assert "<table" not in report


def test_call_command_writes_tsv_summary_for_fuzzy_itd(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(
        ">FLT3 exon 14-15 assay\nAAACCCGGGTTT\n",
        encoding="utf-8",
    )
    report_path = tmp_path / "calls.tsv"

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAAAAACCCGGGCCCGGATTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "AAAAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            f"AAA{reverse_complement('AAACCCGGGCCCGGATTT')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "AAAAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
                "--max-copy-mismatch-rate",
                str(1 / 6),
                "--output-tsv",
                str(report_path),
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""
    rows = list(csv.reader(report_path.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    assert rows[0] == [
        "Call",
        "Status",
        "Filter Reasons",
        "Mode",
        "Max Copy Mismatch Rate",
        "Insertion After Reference Base",
        "Copied Segment Start",
        "Copied Segment End",
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
    assert rows[1][1] == "FAIL"
    assert rows[1][2] == (
        "LOW_MUTANT_FRAGMENT_COUNT;LOW_INFORMATIVE_FRAGMENT_COUNT"
    )
    assert rows[1][3] == "fuzzy"
    assert rows[1][4] == "0.166667"
    assert rows[1][5] == "9"
    assert rows[1][6] == "4"
    assert rows[1][7] == "9"
    assert rows[1][8] == "CCCGGG"
    assert rows[1][11] == "CCCGGA"
    assert rows[1][14:20] == [
        "1",
        "2",
        "0.500000",
        "1",
        "2",
        "0.500000",
    ]
    assert rows[1][20:23] == [
        "2",
        "0.500000",
        "1/2 informative fragments",
    ]
    assert rows[1][32:36] == ["30", "30.000000", "15", "3"]
    assert rows[1][36:41] == [
        "sample",
        "complete",
        "fail",
        "indeterminate",
        "LOW_USABLE_FRAGMENT_COUNT;LOW_MEDIAN_INTERBASE_COVERAGE",
    ]
    assert rows[1][62:64] == ["0", "1"]
    assert rows[1][69:75] == ["1", "0", "0", "0", "1", "0"]
    assert rows[1][75:77] == ["FLT3 exon 14-15 assay", "12"]
    assert rows[1][77] == hashlib.sha256(b"AAACCCGGGTTT").hexdigest()
    assert rows[1][78] == "immediately before insertion"
    assert rows[1][79:82] == ["6", "6", "Yes"]
    assert rows[1][82:87] == [
        "No",
        "0.125000",
        "1.000000",
        "0.050000",
        "3",
    ]
    assert rows[1][87:] == ["0", "0", "."]


def test_adequate_no_call_sample_is_reported_as_qc_passing_negative(
    tmp_path,
    capsys,
) -> None:
    reference = "AAACCCGGGTTT"
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(f">FLT3\n{reference}\n", encoding="utf-8")
    r1_path = tmp_path / "negative_R1.fastq"
    r2_path = tmp_path / "negative_R2.fastq"
    records = "".join(
        f"@fragment-{index}/1\nAAA{reference}\n+\n{'I' * (len(reference) + 3)}\n"
        for index in range(10)
    )
    reverse_records = "".join(
        f"@fragment-{index}/2\nAAA{reference}\n+\n{'I' * (len(reference) + 3)}\n"
        for index in range(10)
    )
    r1_path.write_text(records, encoding="utf-8")
    r2_path.write_text(reverse_records, encoding="utf-8")
    html_path = tmp_path / "negative.html"
    tsv_path = tmp_path / "negative.tsv"

    assert cli.main(
        [
            "--reference",
            str(reference_path),
            "--r1",
            str(r1_path),
            "--r2",
            str(r2_path),
            "--forward-primer",
            "AAA",
            "--reverse-primer",
            "AAA",
            "--min-read-length",
            "12",
            "--output",
            str(html_path),
            "--output-tsv",
            str(tsv_path),
        ]
    ) == 0

    assert capsys.readouterr().out == ""
    report = html_path.read_text(encoding="utf-8")
    assert '<div class="sample-name">ITDiscover Report</div>' in report
    assert '<h2 class="outcome">' not in report
    assert "QC Status:</span><span class=\"status-value status-value--success\">PASS" in report
    assert "No passing ITD was detected." in report
    rows = list(
        csv.reader(tsv_path.read_text(encoding="utf-8").splitlines(), delimiter="\t")
    )
    assert len(rows) == 2
    assert rows[1][0] == "."
    assert rows[1][36:41] == [
        "negative",
        "complete",
        "pass",
        "no passing ITD detected",
        ".",
    ]


def test_cli_can_report_short_out_of_frame_tandem_when_explicitly_enabled(
    tmp_path,
) -> None:
    reference = "GGGATGCCCTACTTT"
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(f">FLT3\n{reference}\n", encoding="utf-8")
    r1_path = tmp_path / "short_R1.fastq"
    r2_path = tmp_path / "short_R2.fastq"
    mutant = "GGGATGCCCACCCTACTTT"
    r1_path.write_text(
        f"@fragment/1\n{mutant}\n+\n{'I' * len(mutant)}\n",
        encoding="utf-8",
    )
    r2_path.write_text(
        f"@fragment/2\n{reverse_complement(reference)}\n+\n"
        f"{'I' * len(reference)}\n",
        encoding="utf-8",
    )
    report_path = tmp_path / "short.tsv"

    assert cli.main(
        [
            "--reference",
            str(reference_path),
            "--r1",
            str(r1_path),
            "--r2",
            str(r2_path),
            "--forward-primer",
            "GGG",
            "--reverse-primer",
            "AAA",
            "--min-read-length",
            "12",
            "--min-insert-length",
            "4",
            "--min-copied-segment-length",
            "3",
            "--no-require-in-frame",
            "--output-tsv",
            str(report_path),
        ]
    ) == 0

    rows = list(
        csv.reader(
            report_path.read_text(encoding="utf-8").splitlines(),
            delimiter="\t",
        )
    )
    assert rows[1][8] == "CCC"
    assert rows[1][11] == "CCCA"
    assert rows[1][79:82] == ["4", "3", "No"]


def test_analysis_error_report_is_indeterminate(tmp_path) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")
    r1_path = tmp_path / "broken_R1.fastq"
    r2_path = tmp_path / "broken_R2.fastq"
    r1_path.write_text("@broken/1\nAAACCC\n+\n", encoding="utf-8")
    r2_path.write_text("", encoding="utf-8")
    report_path = tmp_path / "error.html"

    with pytest.raises(ValueError, match="incomplete"):
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--output",
                str(report_path),
            ]
        )

    report = report_path.read_text(encoding="utf-8")
    assert "Analysis Status" in report
    assert "ERROR" in report
    assert '<h2 class="outcome">' not in report
    assert "ANALYSIS_ERROR" in report
    assert "ITD calling did not complete" in report


def test_call_command_trims_configured_primers(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "GGGTTTAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "GGGTTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            f"{reverse_complement('AAACCCGGGCCCGGGTTTCGTAAA')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            f"{reverse_complement('AAACCCGGGTTTCGTAAA')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "TTT",
                "--reverse-primer",
                "TTTACG",
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""


def test_call_command_silently_filters_reads(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
            "@low-length/1\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
            "@low-quality/1\n"
            "AAACCCGGGTTTAA\n"
            "+\n"
            "!!!!!!!!!!!!!!\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
            "@low-length/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
            "@low-quality/2\n"
            "AAACCCGGGTTTAA\n"
            "+\n"
            "!!!!!!!!!!!!!!\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--min-read-length",
                "13",
                "--min-mean-quality",
                "30",
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""


def test_call_command_rejects_multi_sequence_reference(tmp_path) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">ref1\nAAAA\n>ref2\nCCCC\n", encoding="utf-8")

    with pytest.raises(ValueError, match="exactly one sequence"):
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                "unused_R1.fastq",
                "--r2",
                "unused_R2.fastq",
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
            ]
        )


def test_call_command_writes_concise_html_report(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nTTTAAACCCGGGTTT\n", encoding="utf-8")
    report_path = tmp_path / "reports" / "unique-support.html"

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment-1/1\n"
            "TTTTAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIIII\n"
            "@itd-fragment-2/1\n"
            "TTCAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "TTTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment-1/2\n"
            f"AAA{reverse_complement('TTTAAACCCGGGCCCGGGTTT')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIIIIII\n"
            "@itd-fragment-2/2\n"
            f"AAA{reverse_complement('TTCAAACCCGGGCCCGGGTTT')}\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "AAATTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "T",
                "--reverse-primer",
                "AAA",
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
                "--output",
                str(report_path),
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""
    report = report_path.read_text(encoding="utf-8")
    assert "<title>ITDiscover Report</title>" in report
    assert "<h1>ITDiscover Report</h1>" not in report
    assert '<h2 class="outcome">' not in report
    assert "No passing ITD result is available" in report
    assert "Representative alignment" not in report
    assert "Inserted sequence pileup" not in report
    assert "Reference Sequence SHA-256" not in report
    assert "Concordant Fragments" not in report
    assert "R1 Evidence" not in report
    assert "CALL THRESHOLDS" not in report


def test_html_report_orders_itds_by_support_count_descending(tmp_path) -> None:
    report_path = tmp_path / "ordered-report.html"

    higher_itd = ITD(
        insertion=Insertion(
            read_id="higher-read",
            fragment_id="higher-fragment",
            start=2,
            sequence="CCCGGG",
            direction="forward",
        ),
        copied_segment_start=3,
        copied_segment_sequence="CCCGGG",
    )
    lower_itd = ITD(
        insertion=Insertion(
            read_id="lower-read",
            fragment_id="lower-fragment",
            start=8,
            sequence="TTT",
            direction="forward",
        ),
        copied_segment_start=9,
        copied_segment_sequence="TTT",
    )

    calls = [
        ITDCall(
            itd=lower_itd,
            canonical_allele=CanonicalInsertionAllele(
                start=8,
                sequence="TTT",
            ),
            mutant_fragment_count=1,
            wild_type_fragment_count=9,
        ),
        ITDCall(
            itd=higher_itd,
            canonical_allele=CanonicalInsertionAllele(
                start=2,
                sequence="CCCGGG",
            ),
            mutant_fragment_count=5,
            wild_type_fragment_count=5,
            consolidated_members=(
                ConsolidatedAlleleMember(
                    allele=CanonicalInsertionAllele(
                        start=2,
                        sequence="CCCCGG",
                    ),
                    fragment_count=1,
                    allele_mismatches=1,
                    allele_mismatch_rate=1 / 6,
                    breakpoint_shift=0,
                    breakpoint_shift_rate=0,
                    reason="same-breakpoint sequence error",
                ),
            ),
        ),
    ]
    cli._write_html_report(
        report_path,
        calls,
        [],
    )

    report = report_path.read_text(encoding="utf-8")
    assert report.index("CCCGGG") < report.index("TTT")
    assert "ITD 1" not in report
    assert "ITD 2" not in report
    assert "FLT3 internal tandem duplication" not in report
    assert "Minor-variant consolidation" not in report
    assert "sequence=CCCCGG" not in report
    assert "Inserted sequence pileup" not in report


def test_representative_alignment_prefers_common_exact_support() -> None:
    itd = ITD(
        insertion=Insertion(
            read_id="exact-high",
            fragment_id="fragment",
            start=2,
            sequence="CCCGGG",
            direction="forward",
        ),
        copied_segment_start=3,
        copied_segment_sequence="CCCGGG",
    )
    allele = CanonicalInsertionAllele(start=2, sequence="CCCGGG")

    def representative(
        read_id: str,
        *,
        support_count: int,
        exact_support_count: int,
        mismatches: int,
    ) -> UniqueSupportRepresentative:
        return UniqueSupportRepresentative(
            itd=itd,
            signature=read_id,
            alignment=Alignment(
                read_id=read_id,
                fragment_id=read_id,
                read_sequence="AAACCCGGGCCCGGGTTT",
                aligned_reference="AAACCCGGG------TTT",
                aligned_read="AAACCCGGGCCCGGGTTT",
                direction="forward",
            ),
            canonical_allele=allele,
            support_count=support_count,
            exact_support_count=exact_support_count,
            mismatches=mismatches,
        )

    fuzzy_high = representative(
        "fuzzy-high",
        support_count=10,
        exact_support_count=0,
        mismatches=1,
    )
    exact_low = representative(
        "exact-low",
        support_count=1,
        exact_support_count=1,
        mismatches=0,
    )
    exact_high = representative(
        "exact-high",
        support_count=3,
        exact_support_count=3,
        mismatches=0,
    )

    assert cli._best_representative(
        [fuzzy_high, exact_low, exact_high]
    ) == exact_high


def test_alignment_difference_classes_do_not_color_indels_yellow() -> None:
    itd = ITD(
        insertion=Insertion(
            read_id="read",
            fragment_id="fragment",
            start=8,
            sequence="GGG",
            direction="forward",
        ),
        copied_segment_start=9,
        copied_segment_sequence="GGG",
    )
    deletion = Alignment(
        read_id="deletion-read",
        fragment_id="fragment",
        read_sequence="AAA",
        aligned_reference="AAAC",
        aligned_read="AAA-",
        direction="forward",
    )
    insertion = Alignment(
        read_id="insertion-read",
        fragment_id="fragment",
        read_sequence="AAATC",
        aligned_reference="AAA-C",
        aligned_read="AAATC",
        direction="forward",
    )
    substitution = Alignment(
        read_id="substitution-read",
        fragment_id="fragment",
        read_sequence="AAAT",
        aligned_reference="AAAC",
        aligned_read="AAAT",
        direction="forward",
    )

    assert "diff" not in cli._alignment_difference_classes(deletion, itd)
    assert "diff" not in cli._alignment_difference_classes(insertion, itd)
    assert cli._alignment_difference_classes(substitution, itd)[3] == "diff"


def test_alignment_difference_classes_color_spacers_and_inserted_sequence() -> None:
    itd = ITD(
        insertion=Insertion(
            read_id="read",
            fragment_id="fragment",
            start=2,
            sequence="NNNCCCGGGNN",
            direction="forward",
        ),
        copied_segment_start=3,
        copied_segment_sequence="CCCGGG",
        spacer_prefix="NNN",
        spacer_suffix="NN",
    )
    alignment = Alignment(
        read_id="read",
        fragment_id="fragment",
        read_sequence="AAANNNCCCGGGNNCCCGGGTTT",
        aligned_reference="AAA-----------CCCGGGTTT",
        aligned_read="AAANNNCCCGGGNNCCCGGGTTT",
        direction="forward",
    )

    classes = cli._alignment_difference_classes(alignment, itd)

    assert classes[3:14] == [
        "spacer-region",
        "spacer-region",
        "spacer-region",
        "inserted-region",
        "inserted-region",
        "inserted-region",
        "inserted-region",
        "inserted-region",
        "inserted-region",
        "spacer-region",
        "spacer-region",
    ]


def test_representative_alignment_highlights_all_evidence_categories() -> None:
    itd = ITD(
        insertion=Insertion(
            read_id="read",
            fragment_id="fragment",
            start=2,
            sequence="NNNCCCGGGNN",
            direction="forward",
        ),
        copied_segment_start=3,
        copied_segment_sequence="CCCGGG",
        spacer_prefix="NNN",
        spacer_suffix="NN",
    )
    representative = UniqueSupportRepresentative(
        itd=itd,
        signature="support",
        alignment=Alignment(
            read_id="read",
            fragment_id="fragment",
            read_sequence="AAANNNCCCGGGNNCCCGGGTTA",
            aligned_reference="AAA-----------CCCGGGTTT",
            aligned_read="AAANNNCCCGGGNNCCCGGGTTA",
            direction="forward",
        ),
        canonical_allele=CanonicalInsertionAllele(start=2, sequence="NNNCCCGGGNN"),
        support_count=1,
        exact_support_count=1,
    )

    alignment_html = cli._render_html_representative_alignment(representative)

    assert 'class="tandem-region"' in alignment_html
    assert 'class="inserted-region"' in alignment_html
    assert 'class="spacer-region"' in alignment_html
    assert 'class="diff"' in alignment_html


def test_reference_position_markers_are_1_based_and_ignore_alignment_gaps() -> None:
    markers = cli._render_reference_position_markers("AAA---CCCGGGTTT")

    assert ">10</span>" in markers
    assert 'left: calc(11ch + 12ch)' in markers


def test_call_command_rejects_non_html_output_path(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                "unused_R1.fastq",
                "--r2",
                "unused_R2.fastq",
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--output",
                "report.txt",
            ]
        )

    assert exc_info.value.code == 2
    assert "must end with .html" in capsys.readouterr().err


def test_call_command_rejects_invalid_max_copy_mismatch_rate(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "--reference",
                "reference.fasta",
                "--r1",
                "sample_R1.fastq",
                "--r2",
                "sample_R2.fastq",
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--max-copy-mismatch-rate",
                "-1",
            ]
        )

    assert exc_info.value.code == 2
    assert "value must be between 0 and 1" in capsys.readouterr().err


def test_call_command_filters_without_stdout_report(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--forward-primer",
                "AAA",
                "--reverse-primer",
                "AAA",
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
                "--min-mutant-fragment-count",
                "2",
                "--min-mutant-fragment-fraction",
                "0.6",
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""


def test_alignment_comparison_classes_ignore_leading_gap_shift() -> None:
    baseline = Alignment(
        read_id="baseline",
        fragment_id="baseline",
        read_sequence="GCAATTTAGGT",
        aligned_reference="GCAATTTAGGT",
        aligned_read="GCAATTTAGGT",
        direction="forward",
    )
    shifted = Alignment(
        read_id="shifted",
        fragment_id="shifted",
        read_sequence="AGCAATTTAGGT",
        aligned_reference="-GCAATTTAGGT",
        aligned_read="AGCAATTTAGGT",
        direction="forward",
    )

    classes = cli._alignment_comparison_classes(shifted, baseline)

    assert classes[0] == "insert"
    assert all(css_class is None for css_class in classes[1:])
