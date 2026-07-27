# Alignment walkthrough

This small, full-length illustrative example shows how ITDiscover moves from
read alignments to candidate-specific fragment evidence. The sequences are
deliberately short so every base fits on the page; they demonstrate the same
logic used for full amplicon reads.

## Candidate allele

The reference sequence is:

```text
             1         2
REF ACGTTGCACTGAACTGCTACGATCGTAC
```

The candidate inserts `GAA` after reference base 10. `GAA` is also reference
bases 11–13, immediately after the breakpoint, so this is an ITD candidate.

```text
             1            2
REF ACGTTGCACT---GAACTGCTACGATCGTAC
ALT ACGTTGCACTGAAGAACTGCTACGATCGTAC
              ^^^ inserted GAA
```

For clarity, assume all displayed junction-anchor bases and inserted bases pass
the configured quality filters. In a real run, those quality checks happen
before the evidence state is assigned.

## 1. Opportunities are mutant or wild-type junction observations

An R1 or R2 opportunity is a read that covers the candidate junction well
enough to be classified as either mutant or wild type. The following two
alignments show one of each outcome.

A full-length R1 spans the candidate junction:

```text
             1            2
REF ACGTTGCACT---GAACTGCTACGATCGTAC
R1  ACGTTGCACTGAAGAACTGCTACGATCGTAC
```

The gap in the reference alignment is the `GAA` insertion. R1 therefore
supports the candidate ALT and contributes one R1 mutant observation and one
R1 opportunity.

Its reference-oriented R2 begins downstream of the left junction anchor:

```text
             2
REF GAACTGCTACGATCGTAC
R2  GAACTGCTACGATCGTAC
```

R2 is a valid alignment, but it cannot classify the junction because it does
not cover the bases on both sides of it. It contributes neither an R2 mutant
observation nor an R2 opportunity. The paired fragment is `single-mate` mutant
support, not a conflicting fragment.

A full-length R1 can instead match the reference across the same site:

```text
             1         2
REF ACGTTGCACTGAACTGCTACGATCGTAC
R1  ACGTTGCACTGAACTGCTACGATCGTAC
```

R1 supports the reference junction. It adds one R1 opportunity, but no R1
mutant observation. If no mate supports a conflicting candidate, the paired
fragment is `wild type` evidence for this candidate.

At this point the R1 rate is 1 mutant / 2 opportunities. The fragment-level
counts are 1 mutant and 1 wild type, giving an observed mutant-fragment
fraction of `1 / (1 + 1) = 50%`.

## 2. Concordant mates still count as one fragment

Two full-length mates can both show the candidate insertion:

```text
             1            2
REF ACGTTGCACT---GAACTGCTACGATCGTAC
R1  ACGTTGCACTGAAGAACTGCTACGATCGTAC
R2  ACGTTGCACTGAAGAACTGCTACGATCGTAC
```

Both directions now have a mutant observation and opportunity, but the pair
adds only **one** mutant fragment. It is reported as `concordant` support. The
fragment-level count prevents paired reads from being treated as independent
observations.

## 3. Contradictory mates are excluded from the fraction

An R1 can support the candidate while its R2 fully supports the reference
junction:

```text
             1            2
REF ACGTTGCACT---GAACTGCTACGATCGTAC
R1  ACGTTGCACTGAAGAACTGCTACGATCGTAC

             1         2
REF ACGTTGCACTGAACTGCTACGATCGTAC
R2  ACGTTGCACTGAACTGCTACGATCGTAC
```

The mates make incompatible claims for the same candidate. Their paired
fragment is `conflicting`, not mutant and not wild type. It is excluded from
the observed mutant-fragment-fraction denominator, while remaining visible in
the TSV audit columns.

## 4. Read geometry can leave one direction without opportunities

A read direction that does not reach the junction is not counted as wild type
merely because it lacks the insertion. It cannot distinguish ALT from
reference. This is the reason a report can legitimately show many opportunities
in one direction and zero in the other for an insertion near a read edge.

Direction bias is not assessed until both directions meet the configured
minimum opportunity count. This avoids rejecting an ITD simply because the
assay's read geometry leaves one mate unable to observe that junction.

## 5. An ITD can include a spacer

The copied reference segment does not have to occupy the entire insertion. In
this example, the sequence immediately after reference base 10 starts with
`GAACTG`. The observed insertion is `TGAACTG`:

```text
             1                2
REF ACGTTGCACT-------GAACTGCTACGATCGTAC
R1  ACGTTGCACTTGAACTGGAACTGCTACGATCGTAC

Annotated insertion: [T][GAACTG]
                      spacer prefix | copied reference bases 11–16
```

The seven dashes in the reference alignment correspond to the seven inserted
bases `TGAACTG` in R1. The brackets annotate portions of that sequence; they
are not bases. Because `GAACTG` is an exact adjacent copy, the read is an ITD
observation. The leading `T` is reported as `Spacer Prefix`; it is retained in
the insertion sequence, not discarded or treated as a mismatch.

The extra bases can instead follow the copied segment. For an observed
insertion of `GAACTGT`, the same six-base copied segment is reported with `T`
as `Spacer Suffix`:

```text
             1                2
REF ACGTTGCACT-------GAACTGCTACGATCGTAC
R1  ACGTTGCACTGAACTGTGAACTGCTACGATCGTAC

Annotated insertion: [GAACTG][T]
                      copied reference bases 11–16 | spacer suffix
```

Here, the final `T` is the `Spacer Suffix`.

The copied reference segment must be directly adjacent to the insertion
breakpoint.

## 6. Fuzzy matching tolerates a limited copy difference

Fuzzy matching is an ITD-detection setting, controlled by
`--max-copy-mismatch-rate`. It compares the observed copied segment with the
adjacent reference segment. It is separate from minor-variant consolidation.

Suppose the adjacent reference segment is `GAACTG`, while the observed
insertion is `GAATTG`—one base differs:

```text
             1               2
REF ACGTTGCACT------GAACTGCTACGATCGTAC
R1  ACGTTGCACTGAATTGGAACTGCTACGATCGTAC

Copy (bases 11–16) GAACTG
Observed insertion GAATTG
                   ||| ||
                      ^ one C→T difference
```

The mismatch rate is `1 / 6 = 0.1667`. Exact mode (the default,
`--max-copy-mismatch-rate 0`) rejects this as an ITD. A fuzzy threshold of at
least `0.1667` accepts it, provided the read and local evidence also pass their
quality filters. The report retains the observed insertion sequence and the
reference copied-segment sequence, so the difference remains inspectable.

Fuzzy matching does not turn an arbitrary insertion into an ITD: the matching
reference sequence must be adjacent to the insertion breakpoint, and the
matched segment must meet `--min-copied-segment-length`.

## 7. Similar weak variants can be consolidated into one call

Consolidation is optional (`--consolidate-minor-itd-variants`) and happens only
after ITD alleles have been detected. It is an error-suppression heuristic for
weak, compatible minor observations; it is not enabled by default because it
changes reported allele identity and does not prove biological equivalence.

This is different from canonical normalization. Equivalent alignment-gap
placements of the same ALT are normalized to one allele before consolidation.
Consolidation instead considers non-equivalent detected alleles, which can be
technical artefacts but can also be genuine distinct biological variants.

Consider two full-length, already-detected alleles at the same breakpoint:

```text
             2
REF GAACTGCTACGA
ALT GAACTGCTACGA   anchor: 100 supporting fragments
ALT GAATTGCTACGA   minor: 4 supporting fragments
       ^ one difference
```

The two 12-base ALT sequences differ at one position, an allele mismatch rate
of `1 / 12 = 0.0833`. With the default consolidation thresholds, the minor
allele can be assigned to the anchor because all of these conditions hold:

| Check | Example | Default limit |
|---|---:|---:|
| Anchor support | 100 fragments | at least 3 |
| Minor / anchor support | `4 / 100 = 0.04` | at most 0.05 |
| Allele mismatch rate | `1 / 12 = 0.0833` | at most 0.125 |
| Breakpoint shift rate | `0 / 12 = 0` | at most 1 |

Because this minor allele has a copy difference, it must first be detected in
fuzzy mode with a copy-mismatch threshold of at least `1 / 12`. Consolidation
does not discover it in exact mode; it only decides whether to assign an
already-detected minor allele to an anchor.

Consolidation can also assign a weak allele with a nearby, shifted insertion
breakpoint. For a 12-base insertion, a one-base shift has a breakpoint-shift
rate of `1 / 12 = 0.0833`, which is within the default limit of 1. The shifted
allele must still be fully observed, have the same insertion length as the
anchor, satisfy the full alternate-sequence mismatch and support-ratio limits,
and have one uniquely best compatible anchor. A nearby shifted breakpoint may
be a genuine distinct allele or subclone, so enable this option only when this
error-suppression trade-off is appropriate for the assay. Such an assignment is
reported as a nearby-breakpoint local-haplotype match in the TSV.

Spacers are allowed. The required equal length is the total observed insertion
length, including any spacer bases. Consolidation compares the complete
insertion allele and does not require the minor and anchor to have the same
copied-segment or spacer annotation.

The final anchor call then receives the four minor observations in addition to
its own support. The TSV retains the raw minor support and the reason for the
assignment in the consolidated-minor columns, rather than silently hiding the
minor allele.

For contrast, a minor allele supported by 6 fragments would not meet the
default support ratio (`6 / 100 = 0.06`), so it would remain a separate call.
Likewise, an allele without a uniquely best compatible anchor is not assigned.
Assignments are direct; a minor allele cannot become an anchor for another
minor allele.

## From these examples to a passing call

These examples use too few fragments to pass the default minimums. In a normal
run, ITDiscover aggregates the candidate-specific evidence states across all
fragments, then applies the configured mutant-count, informative-depth,
mutant-fraction, ambiguity, and—when evaluable—directional filters. See
[Algorithm and evidence model](algorithm.md) and
[Inputs, outputs, QC, and coordinate conventions](input-output.md) for the
remaining details.
