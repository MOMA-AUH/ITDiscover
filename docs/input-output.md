# Inputs, outputs, QC, and coordinate conventions

## Inputs

ITDiscover requires one reference amplicon FASTA sequence and synchronized R1
and R2 FASTQ files. Paired records must have matching fragment identifiers.

Primer trimming is required. Copy `--forward-primer` directly from the start
of raw R1 reads and `--reverse-primer` directly from the start of raw R2 reads.
Do **not** reverse-complement the R2 primer before supplying it. ITDiscover
reverse-complements R2 to the forward-reference orientation, then trims the
reverse complement of the supplied R2 primer from the oriented read's 3′ end.

Reads are screened for trimmed length, mean quality, reference identity, and
on-target fraction. Candidate evidence is additionally screened using the
junction-anchor and inserted-base quality settings. Use `itdiscover --help` to
view every threshold and its default.

## Reports

`--output` writes a concise HTML summary for review. `--output-tsv` writes a
tab-separated audit report with passing calls and filtered candidates, all
settings, QC metrics, reference identity, and fragment-evidence counts.

The report separates these concepts:

| Field | Values | Meaning |
|---|---|---|
| Analysis status | `complete`, `error` | Whether analysis completed. |
| QC status | `pass`, `warn`, `fail` | Whether sample-level QC passed. |
| Outcome | `ITD detected`, `no passing ITD detected`, `indeterminate` | Final sample interpretation. |

An empty passing-call list is not automatically negative. A sample becomes `no
passing ITD detected` only after sample QC passes. A sample with filtered
candidates but no passing call can be reported as a QC warning or
indeterminate, depending on the evidence and QC result.

## Default sample QC

The default QC criteria require:

- at least 10 usable fragments after preprocessing;
- at least one alignment-passing read in each direction;
- an alignment pass fraction of at least 80%;
- median inter-base fragment coverage of at least 10; and
- at least 80% primer retention in each direction when a primer is configured.

These are configurable research defaults. They do not replace assay-specific
validation.

## Coordinates

All report coordinates are reference-local and 1-based. An insertion coordinate
names the reference base immediately before the insertion; an insertion before
the first base is therefore reported as after base `0`.

Copied-segment start and end coordinates are 1-based and inclusive. A copied
segment can be immediately before or immediately after the insertion in this
coordinate representation. The TSV also records the complete FASTA header,
reference length, and SHA-256 digest so that coordinates remain tied to a
specific reference sequence.

## Interpretation boundary

The observed mutant-fragment fraction is a post-filter evidence measure. It
does not deduplicate PCR duplicates or UMIs, correct amplification bias, or
estimate a clinically validated VAF or allelic ratio. Do not treat it as
interchangeable with a VAF or capillary-electrophoresis allelic ratio from a
different library method.

For the evidence definitions behind the report columns, see the
[algorithm and evidence model](algorithm.md). For a worked set of alignments,
see the [alignment walkthrough](alignment-walkthrough.md).
