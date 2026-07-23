import csv
import hashlib
from pathlib import Path

import pytest

import itdiscover
import itdiscover.cli as cli
from itdiscover.alleles import CanonicalInsertionAllele
from itdiscover.calls import ITDCall, InsertSequenceSupport, UniqueSupportRepresentative
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
    args = cli.build_parser().parse_args(
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
        ]
    )

    assert args.max_directional_mutant_fraction_share == 0.8
    assert args.min_directional_opportunities == 7
    assert args.min_copied_segment_length == 9


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
    assert passing["QC Status"] == "warn"
    assert passing["Insertion Sequence"] == "AGAGAATATGAATAT"
    insertion_coordinate = (
        "Insertion After Reference Base (0-based; -1=before first)"
    )
    assert passing[insertion_coordinate] == "78"
    assert passing["Copied Segment Start (0-based)"] == "79"
    assert passing["Copied Segment End (0-based, inclusive)"] == "93"
    assert passing["Copied Segment Location"] == "immediately after insertion"
    assert passing["Mutant Fragment Count"] == "3"
    assert passing["Wild-type Fragment Count"] == "9"
    assert passing["Informative Fragment Count"] == "12"
    assert passing["Observed Mutant-fragment Fraction"] == "0.250000"
    assert filtered["Status"] == "FAIL"
    assert filtered["Filter Reasons"] == "LOW_MUTANT_FRAGMENT_COUNT"


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
                "--max-mismatches",
                "1",
                "--output",
                str(report_path),
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == ""
    report = report_path.read_text(encoding="utf-8")
    assert "<title>ITDiscover Report</title>" in report
    assert "<h1>ITDiscover Report</h1>" in report
    assert "CALL THRESHOLDS" in report
    assert "Min mutant fragments" in report
    assert "Min informative fragments" in report
    assert "Min mutant-fragment fraction" in report
    assert "Max directional mutant-fraction share" in report
    assert "Min opportunities per direction" in report
    assert "Min insert length" in report
    assert "Min copied-segment length" in report
    assert "Require in-frame insertions" in report
    assert "VAF" not in report
    assert "Max mismatches" in report
    assert ">1<" in report
    assert "Representative alignment" in report
    assert "Concordant Fragments" in report
    assert "Single-mate Fragments" in report
    assert "Conflicting Fragments" in report
    assert "Unresolved Fragments" in report
    assert "Wild-type Fragments" in report
    assert "Not-informative Fragments" in report
    assert "R1 Evidence" in report
    assert "R2 Evidence" in report
    assert "50.0% (1/2 opportunities)" in report
    assert "copied reference segment" in report
    assert "inserted sequence" in report
    assert "spacer sequence" in report
    assert "mismatches" in report
    assert "sky blue" not in report
    assert "teal green" not in report
    assert "orange" not in report
    assert 'class="legend-chip tandem-region"' in report
    assert 'class="tandem-region"' in report
    assert 'class="inserted-region' in report
    assert "<strong>itd-fragment/1</strong>" in report
    assert '<div class="signature">' not in report
    assert "mismatches 1" not in report
    assert "support pattern count 1" not in report
    assert '<span class="support-meta">fragment' not in report
    assert "Inserted sequence pileup" in report
    assert "<th>Inserted sequence</th><th>Mismatches</th><th>Count</th>" in report
    assert "CCCGGG" in report
    assert 'class="inserted-region insert-mismatch"' in report


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
                "--max-mismatches",
                "1",
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
    ]
    assert rows[1][1] == "FAIL"
    assert rows[1][2] == (
        "LOW_MUTANT_FRAGMENT_COUNT;LOW_INFORMATIVE_FRAGMENT_COUNT"
    )
    assert rows[1][3] == "fuzzy"
    assert rows[1][4] == "1"
    assert rows[1][5] == "8"
    assert rows[1][6] == "3"
    assert rows[1][7] == "8"
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
    assert rows[1][34:39] == [
        "sample",
        "complete",
        "fail",
        "indeterminate",
        "LOW_USABLE_FRAGMENT_COUNT;LOW_MEDIAN_INTERBASE_COVERAGE",
    ]
    assert rows[1][60:62] == ["0", "1"]
    assert rows[1][67:73] == ["1", "0", "0", "0", "1", "0"]
    assert rows[1][73:75] == ["FLT3 exon 14-15 assay", "12"]
    assert rows[1][75] == hashlib.sha256(b"AAACCCGGGTTT").hexdigest()
    assert rows[1][76] == cli.COORDINATE_CONVENTION
    assert rows[1][77] == "immediately before insertion"
    assert rows[1][-3:] == ["6", "6", "Yes"]


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
    assert "Sample Result and QC" in report
    assert "negative" in report
    assert "no passing ITD detected" in report
    assert "No ITD candidates were called." in report
    rows = list(
        csv.reader(tsv_path.read_text(encoding="utf-8").splitlines(), delimiter="\t")
    )
    assert len(rows) == 2
    assert rows[1][0] == "."
    assert rows[1][34:39] == [
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
    assert rows[1][-3:] == ["4", "3", "No"]


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
    assert "error" in report
    assert "indeterminate" in report
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


def test_call_command_writes_unique_support_alignment_html_report(tmp_path, capsys) -> None:
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
    assert "<h1>ITDiscover Report</h1>" in report
    assert 'class="legend-chip diff"' in report
    assert "copied reference segment" in report
    assert "mismatches" in report
    assert "sky blue" not in report
    assert "teal green" not in report
    assert "orange" not in report
    assert "<h2>ITD 1</h2>" not in report
    assert "<h2>ITD 2</h2>" not in report
    assert "Insertion After Reference Base (0-based)" in report
    assert "Copied Segment Start (0-based)" in report
    assert "Copied Segment End (0-based, inclusive)" in report
    assert "Copied Segment Location" in report
    assert "Immediately after the insertion" in report
    assert "Tandem Orientation" not in report
    assert "Reference and Coordinates" in report
    assert "Reference FASTA Header" in report
    assert "Reference Sequence SHA-256" in report
    assert cli.COORDINATE_CONVENTION in report
    assert "Mutant Fragments" in report
    assert "Observed Mutant-fragment Fraction" in report
    assert "66.7% (2/3 informative fragments)" in report
    assert "VAF" not in report
    assert "support pattern count 1" not in report
    assert "mismatches 0" not in report
    assert '<div class="signature">' not in report
    assert "Inserted sequence pileup" in report
    assert "CCCGGG" in report
    assert '<span class="diff">C</span>' in report


def test_unique_support_report_orders_itds_by_support_count_descending(tmp_path) -> None:
    report_path = tmp_path / "ordered-report.html"

    higher_alignment = Alignment(
        read_id="higher-read",
        fragment_id="higher-fragment",
        read_sequence="AAACCCGGGCCCGGGTTT",
        aligned_reference="AAACCCGGG------TTT",
        aligned_read="AAACCCGGGCCCGGGTTT",
        direction="forward",
    )
    lower_alignment = Alignment(
        read_id="lower-read",
        fragment_id="lower-fragment",
        read_sequence="AAACCCGGGTTT",
        aligned_reference="AAACCCGGGTTT",
        aligned_read="AAACCCGGGTTT",
        direction="forward",
    )

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
        ),
    ]
    representatives = [
        UniqueSupportRepresentative(
            itd=lower_itd,
            signature="lower[sig]",
            alignment=lower_alignment,
            canonical_allele=CanonicalInsertionAllele(
                start=8,
                sequence="TTT",
            ),
            support_count=1,
            exact_support_count=1,
            mismatches=0,
            insert_sequence_supports=(
                InsertSequenceSupport(sequence="TTT", support_count=1, mismatches=0),
            ),
        ),
        UniqueSupportRepresentative(
            itd=higher_itd,
            signature="higher[sig]",
            alignment=higher_alignment,
            canonical_allele=CanonicalInsertionAllele(
                start=2,
                sequence="CCCGGG",
            ),
            support_count=5,
            exact_support_count=5,
            mismatches=0,
            insert_sequence_supports=(
                InsertSequenceSupport(
                    sequence="CCCGGG",
                    support_count=5,
                    mismatches=0,
                ),
            ),
        ),
    ]

    cli._write_unique_support_alignment_html_report(report_path, calls, representatives)

    report = report_path.read_text(encoding="utf-8")
    assert report.index("<strong>higher-read</strong>") < report.index(
        "<strong>lower-read</strong>"
    )


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


def test_call_command_rejects_negative_max_mismatches(capsys) -> None:
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
                "--max-mismatches",
                "-1",
            ]
        )

    assert exc_info.value.code == 2
    assert "value must not be negative" in capsys.readouterr().err


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
