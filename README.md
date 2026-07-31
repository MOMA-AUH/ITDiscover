# ITDiscover

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/itdiscover?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/itdiscover) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/itdiscover?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/itdiscover)

ITDiscover discovers internal tandem duplications (ITDs) from paired-end
amplicon sequencing reads.

It accepts a single-reference amplicon FASTA plus paired FASTQ reads, trims the
assay primers, and writes a concise HTML result summary and/or a complete TSV
audit report.

## Install

Install from the `MOMA-AUH` Conda channel:

```bash
conda install MOMA-AUH::itdiscover
```

## Run

Provide the reference, paired FASTQ files, and primer sequences:

```bash
itdiscover \
  --reference reference.fa \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --forward-primer FORWARD_PRIMER \
  --reverse-primer REVERSE_PRIMER \
  --sample-id sample \
  --output report.html \
  --output-tsv report.tsv
```

Copy each primer directly from the beginning (5′ end) of the corresponding raw
FASTQ sequence:

| Option | Supply |
|---|---|
| `--forward-primer` | The sequence at the start of raw R1 reads. |
| `--reverse-primer` | The sequence at the start of raw R2 reads—**not** its reverse complement. |

For example, if raw R2 reads begin `CTTTCAGC...`, supply
`--reverse-primer CTTTCAGC`. ITDiscover reverse-complements R2 to the
forward-reference orientation internally, then removes the reverse complement
of that supplied primer from the oriented read's 3′ end.

Use `itdiscover --help` for all options and defaults.

The HTML report is intended for quick review. The TSV retains calls, filtered
candidates, QC, thresholds, reference identity, and evidence counts.

## Included example

The repository includes a synthetic FLT3 amplicon example:

```bash
mkdir -p example-output

itdiscover \
  --reference tests/data/synthetic_flt3/reference.fa \
  --r1 tests/data/synthetic_flt3/synthetic_R1.fastq \
  --r2 tests/data/synthetic_flt3/synthetic_R2.fastq \
  --forward-primer GCAATTTAGGTATGAAAGCCAGC \
  --reverse-primer CTTTCAGCATTTTGACGGCAACC \
  --sample-id synthetic-flt3 \
  --output example-output/report.html \
  --output-tsv example-output/calls.tsv
```

The primary call has 3 mutant fragments among 12 informative fragments (25%).
This is an observed post-filter fragment fraction, not a validated VAF or
allelic ratio.

## Important limitations

ITDiscover is assay-specific software. Its defaults, thresholds, and
reported fraction require validation for any intended workflow. It does not
deduplicate PCR reads or UMIs, correct amplification bias, or provide a
clinically validated VAF or diagnostic result.

## Documentation

- [Algorithm and evidence model](docs/algorithm.md)
- [Step-by-step alignment walkthrough](docs/alignment-walkthrough.md)
- [Inputs, outputs, QC, and coordinate conventions](docs/input-output.md)

## Development

Run the test suite in the development environment:

```bash
conda env create -f environment.yml
conda run -n itdiscover-dev pytest
```
