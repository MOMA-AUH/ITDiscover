import csv

import pytest

import itdiscover
import itdiscover.cli as cli
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
            "AAACCCGGGCCCGGATTT\n"
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
    assert "Min support count" in report
    assert "Min coverage" in report
    assert "Min VAF" in report
    assert "Max mismatches" in report
    assert ">1<" in report
    assert "Representative alignment" in report
    assert "tandem sequence" in report
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
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")
    report_path = tmp_path / "calls.tsv"

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAACCCGGGCCCGGATTT\n"
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
    assert rows[1][1] == "PASS"
    assert rows[1][3] == "fuzzy"
    assert rows[1][4] == "1"
    assert rows[1][5] == "8"
    assert rows[1][6] == "3"
    assert rows[1][7] == "8"
    assert rows[1][8] == "CCCGGG"
    assert rows[1][11] == "CCCGGA"


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
            "TTTAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@itd-fragment-2/1\n"
            "TTCAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment-1/2\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
            "@itd-fragment-2/2\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "TTTAAACCCGGGTTT\n"
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
    assert "tandem sequence" in report
    assert "mismatches" in report
    assert "sky blue" not in report
    assert "teal green" not in report
    assert "orange" not in report
    assert "<h2>ITD 1</h2>" not in report
    assert "<h2>ITD 2</h2>" not in report
    assert "Insertion Start" in report
    assert "Support Count" in report
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
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
    )
    lower_itd = ITD(
        insertion=Insertion(
            read_id="lower-read",
            fragment_id="lower-fragment",
            start=8,
            sequence="TTT",
            direction="forward",
        ),
        tandem_start=9,
        tandem_sequence="TTT",
        orientation="downstream",
    )

    calls = [
        ITDCall(itd=lower_itd, support_count=1, coverage=10, vaf=0.1),
        ITDCall(itd=higher_itd, support_count=5, coverage=10, vaf=0.5),
    ]
    representatives = [
        UniqueSupportRepresentative(
            itd=lower_itd,
            signature="lower[sig]",
            alignment=lower_alignment,
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
        tandem_start=9,
        tandem_sequence="GGG",
        orientation="downstream",
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
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
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
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
                "--min-support-count",
                "2",
                "--min-vaf",
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
