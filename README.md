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
By default, `PASS` requires at least 3 mutant fragments, 10 informative
fragments, and an observed mutant-fragment fraction of 0.01. These are
conservative research defaults, not clinically validated thresholds.

Evidence is classified separately for each candidate ALT allele and fragment:

- `mutant`: high-quality evidence for the candidate ALT;
- `wild type`: high-quality evidence for the reference junction;
- `conflicting`: the mates support incompatible alleles;
- `unresolved`: the fragment spans the junction but cannot be assigned
  confidently, including candidate evidence rejected by local quality filters;
- `not informative`: the fragment does not cover the candidate junction.

Mutant and wild-type evidence use the same configured junction anchors and
quality threshold. The reported **observed mutant-fragment fraction** is
`mutant / (mutant + wild type)`, displayed with both counts, for example
`12.3% (37/301 informative fragments)`. Conflicting, unresolved, and
not-informative fragments do not enter the fraction but are always reported.
If conflicting plus unresolved fragments outnumber informative fragments, the
candidate fails and a sample without another passing call is indeterminate.

Mate evidence is reconciled before these states are assigned. `concordant`
fragments have compatible candidate evidence from both mates. `single-mate`
fragments have one mutant mate without high-quality contradictory wild-type
evidence from the other mate. Both are subcategories of `mutant` and each
fragment is counted only once.

This post-filter fragment fraction does not collapse PCR duplicates or UMIs,
correct amplification bias, or estimate a clinically validated VAF or allelic
ratio. It must not be interpreted as interchangeable with capillary
electrophoresis allelic ratio or a VAF from another library method.

FASTQ-derived evidence is also screened for read-to-reference identity and
on-target fraction, Phred quality across the inserted sequence and three-base
junction anchors, standard Illumina adapter motifs, and directional read
imbalance.
Use `itdiscover --help` to inspect or override every evidence threshold. Raw
alignment score filtering and rejection of multiply optimal alignments are
available but disabled by default because suitable cutoffs and equivalent-gap
normalization are assay dependent.

## Sample result and QC

Every CLI report separates three concepts:

- `analysis_status`: `complete` or `error`;
- `qc_status`: `pass`, `warn`, or `fail`; and
- `outcome`: `ITD detected`, `no passing ITD detected`, or `indeterminate`.

An empty call list is not automatically negative. `no passing ITD detected` is
reported only after sample QC passes. By default, QC requires at least 10 usable
fragments, at least one alignment-passing read from each direction, an alignment
pass fraction of at least 80%, median inter-base fragment coverage of at least
10, and at least 80% primer retention in each direction for which a primer was
configured. These research defaults are configurable with the `--min-*` QC
options and require assay-specific validation.
An otherwise adequate sample is marked `warn` when filtered ITD candidates are
present, while the outcome continues to state that no *passing* ITD was detected.

HTML and TSV reports retain input fragment/read counts, primer failures and
retention by direction, length and quality failures, preprocessing-passing
reads, usable fragments, alignment retention and direction counts, inter-base
coverage min/median/max, passing call count, filtered candidate count, QC
reasons, thresholds, and analysis errors. A malformed input can therefore be
distinguished from an inadequate assay and from an adequate no-call sample.

Reports identify the reference by its complete FASTA header, sequence length,
and SHA-256 digest. All reported coordinates are zero-based and local to that
reference. The insertion coordinate is the reference base immediately before
the insertion (`-1` means before the first base); tandem start is zero-based and
tandem end is zero-based and inclusive. Tandem orientation is `upstream` when
the copied reference interval ends at the insertion coordinate and `downstream`
when it begins at the following reference base.

## Event-length and frame policy

`--min-insert-length` controls the shortest observed insertion considered.
`--min-tandem-length` independently controls the shortest copied adjacent
reference tract accepted as an ITD. When omitted, the tandem minimum follows
the insertion minimum, so `--min-insert-length 3` can call a genuine 3-base
duplication instead of retaining a hidden 6-base floor.

Fully observed insertions are required to have total inserted length divisible
by three by default. This narrow FLT3-ITD policy is explicit and can be disabled
with `--no-require-in-frame` when exploratory or assay-specific analysis should
retain out-of-frame candidates. Read-edge/partial observations cannot establish
the full event frame and remain handled as partial evidence. The insertion
minimum, tandem minimum, and frame policy are recorded in HTML and TSV output.
Changing these settings changes the reportable candidate space; it does not by
itself validate short or out-of-frame events for clinical interpretation.

Primer trimming is optional and can be enabled with `--forward-primer` and
`--reverse-primer`. Supply each primer in the orientation in which it occurs at
the 5′ end of its raw FASTQ read: the reverse-primer sequence is therefore
reverse-complemented internally before it is trimmed from the 3′ end of the
reference-oriented R2 read.

For an insertion to be called an ITD, its copied reference tract must be immediately
adjacent to the insertion breakpoint. Extra inserted bases may flank the copied tract
and are reported as spacer sequence.
