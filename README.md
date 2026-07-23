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

## Worked FLT3 example

The repository includes synthetic paired-end reads constructed from a real
329-base human FLT3 amplicon reference sequence. The locus is
[human FLT3](https://www.ncbi.nlm.nih.gov/gene/2322), and the 11F/12R primer
pair is a published assay that produces a 329-base wild-type product
([method description](https://www.nature.com/articles/s41467-026-68582-2)).
From a source checkout, run:

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

The expected primary result is:

| Report field | Expected value |
|---|---|
| Analysis status | `complete` |
| QC status | `pass` |
| Outcome | `ITD detected` |
| Inserted sequence | `AGAGAATATGAATAT` |
| Insertion coordinate | after reference base 78 |
| Copied reference segment | bases 79–93, immediately after the insertion |
| Mutant fragments | 3 |
| Wild-type fragments | 9 |
| Informative fragments | 12 |
| Observed mutant-fragment fraction | 25.0% (`3/12`) |

In plain language: an ITD is detected, and 3 of the 12 fragments that can
confidently distinguish this allele from wild type support the ITD. A second
candidate that fails the minimum of three mutant fragments remains visible for
auditability without downgrading the passing call. The 25.0% value is an
observed post-filter fragment fraction, not a clinically validated VAF or
allelic ratio.
The reads are synthetic and the default thresholds are research defaults, so
this example demonstrates interpretation of ITDiscover output rather than a
validated clinical result.

The `--output` flag writes an HTML report with one representative alignment per
called duplication. The `--output-tsv` flag writes the same call and QC facts
as tab-separated data.

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

Weak sequence-error variants can optionally be consolidated into a dominant
ITD allele with `--consolidate-minor-itd-variants`. Consolidation is disabled
by default because it changes reported allele identity. ITD detection occurs
first: `--max-copy-mismatches` controls only how many mismatches are allowed
between the copied segment and the reference. Consolidation is a separate,
advanced second step that compares already-detected alleles. A minor allele is
assigned directly to one uniquely best-supported anchor only when:

- both alleles are fully observed and have the same insertion length;
- their complete ALT sequences differ by no more than
  `--consolidation-max-allele-mismatches` positional mismatches (default 1);
- their breakpoints differ by no more than
  `--consolidation-max-breakpoint-shift` bases (default 6);
- the anchor has at least
  `--consolidation-min-anchor-fragment-count` evidence-passing raw supporting
  fragments
  (default 3); and
- minor support is at most
  `--consolidation-max-minor-support-ratio` of anchor support (default 0.05).

Assignments are never chained through another minor allele. Fragment and mate
consensus is recalculated after assignment, so a fragment is still counted at
most once. HTML and TSV reports record the settings and list every absorbed
allele with its evidence-passing raw fragment support, sequence distance,
breakpoint shift, and reason. This makes aggressive exploratory settings
auditable rather than silently hiding minor observations.

R1/R2 direction bias is evaluated using observation opportunities rather than
raw mutant counts. For each direction, an opportunity is a read that covers the
candidate junction well enough to be classified as high-quality mutant or wild
type. The direction-specific mutant fraction is
`mutant observations / opportunities`. ITDiscover divides the larger of the
two directional fractions by their sum and compares that share with the
configured maximum. It does this only when both directions have at least the
configured number of opportunities. Therefore, an ITD near a read edge is not
called direction-biased merely because only one read direction can reach its
junction.

This post-filter fragment fraction does not collapse PCR duplicates or UMIs,
correct amplification bias, or estimate a clinically validated VAF or allelic
ratio. It must not be interpreted as interchangeable with capillary
electrophoresis allelic ratio or a VAF from another library method.

FASTQ-derived evidence is also screened for read-to-reference identity and
on-target fraction, Phred quality across the inserted sequence and three-base
junction anchors, and directional read imbalance.
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
An otherwise adequate sample with no passing call is marked `warn` when filtered
ITD candidates are present, while the outcome continues to state that no
*passing* ITD was detected. Minor filtered alternatives do not downgrade an
otherwise adequate sample that has a passing ITD call; their count and details
remain available in the reports.

HTML and TSV reports retain input fragment/read counts, primer failures and
retention by direction, length and quality failures, preprocessing-passing
reads, usable fragments, alignment retention and direction counts, inter-base
coverage min/median/max, passing call count, filtered candidate count, QC
reasons, thresholds, and analysis errors. A malformed input can therefore be
distinguished from an inadequate assay and from an adequate no-call sample.

Reports identify the reference by its complete FASTA header, sequence length,
and SHA-256 digest. All reported coordinates are zero-based and local to that
reference. The insertion coordinate is the reference base immediately before
the insertion (`-1` means before the first base). The copied segment start is
zero-based and its end is zero-based and inclusive. Reports state directly
whether that segment is immediately before or immediately after the insertion.
This is a description of the displayed coordinates, not a claim about the
biological direction of copying. In repetitive sequence, the same ALT sequence
can permit more than one equivalent gap placement. ITDiscover groups those
placements as one canonical allele before counting fragments; the
before/after label only explains the selected report representation.

For example, the bundled FLT3 amplicon example contains an insertion after
reference base 78. The inserted sequence `AGAGAATATGAATAT` copies reference
bases 79–93, so the report says that the copied segment is immediately after
the insertion. These coordinates and bases come from
`tests/data/synthetic_flt3/reference.fa`, not from a toy repeat.

## Event-length and frame policy

`--min-insert-length` controls the shortest observed insertion considered.
`--min-copied-segment-length` independently controls the shortest copied
adjacent reference segment accepted as an ITD. When omitted, the copied segment
minimum follows the insertion minimum, so `--min-insert-length 3` can call a
genuine 3-base duplication instead of retaining a hidden 6-base floor.

Fully observed insertions are required to have total inserted length divisible
by three by default. This narrow FLT3-ITD policy is explicit and can be disabled
with `--no-require-in-frame` when exploratory or assay-specific analysis should
retain out-of-frame candidates. Read-edge/partial observations cannot establish
the full event frame and remain handled as partial evidence. The insertion
minimum, copied-segment minimum, and frame policy are recorded in HTML and TSV
output.
Changing these settings changes the reportable candidate space; it does not by
itself validate short or out-of-frame events for clinical interpretation.

Primer trimming is required. Supply `--forward-primer` and `--reverse-primer`
in the orientation in which each occurs at the 5′ end of its raw FASTQ read:
the reverse-primer sequence is therefore reverse-complemented internally before
it is trimmed from the 3′ end of the reference-oriented R2 read.

For an insertion to be called an ITD, its copied reference tract must be immediately
adjacent to the insertion breakpoint. Extra inserted bases may flank the copied tract
and are reported as spacer sequence.
