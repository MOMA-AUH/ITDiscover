# ITDiscover

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/itdiscover?cacheSeconds=300)](https://anaconda.org/MOMA-AUH/itdiscover) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/itdiscover?cacheSeconds=300)](https://anaconda.org/MOMA-AUH/itdiscover)

Tool for discovering FLT3 ITDs from amplicon sequencing of AML samples

## Installation

The recommended way to install **ITDiscover** is via [conda](https://docs.conda.io/), using the `MOMA-AUH` channel:

```bash
conda install MOMA-AUH::itdiscover
```

## Usage

```bash
itdiscover --help
itdiscover --version
```

Example:

```bash
itdiscover \
  --reference reference.fasta \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --forward-primer GGGTTT \
  --reverse-primer AAACCC \
  --output report.html
```

The `--output` flag writes an HTML report with one representative alignment per called duplication.
The `--output-tsv` flag writes a tab-separated summary of the called ITDs.
Primer trimming is optional and can be enabled with `--forward-primer` and
`--reverse-primer`. Supply each primer in the orientation in which it occurs at
the 5′ end of its raw FASTQ read: the reverse-primer sequence is therefore
reverse-complemented internally before it is trimmed from the 3′ end of the
reference-oriented R2 read.

For an insertion to be called an ITD, its copied reference tract must be immediately
adjacent to the insertion breakpoint. Extra inserted bases may flank the copied tract
and are reported as spacer sequence.
