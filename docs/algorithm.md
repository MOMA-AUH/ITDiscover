# Algorithm and evidence model

This document describes the ITDiscover calling model. It is intended to make
the reported evidence counts and filter reasons auditable; it is not a
clinical-validation specification.

## Processing overview

1. Read paired FASTQ records and retain their shared fragment identifier.
2. Orient R1 and R2 to the forward reference, trim the configured primers, and
   apply read-length and mean-quality filters.
3. Align retained reads to the one-sequence reference amplicon and apply
   identity and on-target filters.
4. Extract insertions and apply local junction-anchor and insert-base quality
   requirements.
5. Identify an insertion as an ITD when a copied reference segment is adjacent to
   the insertion breakpoint and satisfies the selected exact or fuzzy copy
   rule.
6. Normalize equivalent gap placements, reconcile mate evidence by fragment,
   then apply call and sample-QC filters.

An optional final consolidation step is an error-suppression heuristic that can
assign a weak, compatible minor ITD allele to a stronger allele. It is disabled
by default because it changes the reported allele identity and does not
establish that the two observations are biologically the same allele.

## Candidate-specific fragment evidence

Evidence is determined separately for every candidate allele. One fragment can
occupy exactly one of these evidence states for that candidate:

| State | Meaning |
|---|---|
| `mutant` | High-quality evidence supports the candidate ALT junction. |
| `wild type` | High-quality evidence supports the reference junction. |
| `conflicting` | The two mates support incompatible alleles, or a candidate and high-quality wild type. |
| `unresolved` | The junction is spanned but cannot be assigned confidently, including locally rejected candidate evidence. |
| `not informative` | Neither mate can classify this candidate junction. |

`mutant` is further reported as `concordant` when compatible candidate evidence
is present in both mates, or `single-mate` when one mate supports the candidate
and the other does not provide high-quality contradictory wild-type evidence.
Each fragment is counted once, even when both mates support it.

## Informative fraction and opportunities

The observed mutant-fragment fraction is:

```text
mutant fragments / (mutant fragments + wild-type fragments)
```

Conflicting, unresolved, and not-informative fragments are reported but do not
enter that denominator.

Directional evidence uses a different denominator. An R1 or R2 opportunity is
a fragment whose read in that direction can classify the candidate junction as
high-quality mutant or wild type:

```text
directional mutant fraction = directional mutant observations / directional opportunities
```

This distinction matters near a read edge: a direction that cannot reach the
junction has zero opportunities, rather than contributing apparent
wild-type support. Direction-bias filtering is performed only when both R1
and R2 have at least `--min-directional-opportunities` opportunities.

## Default calling filters

By default, a candidate needs at least 3 mutant fragments, 10 informative
fragments, and an observed mutant-fragment fraction of 1%. A candidate also
fails when conflicting plus unresolved evidence exceeds mutant plus wild-type
evidence. Directional imbalance is evaluated using the opportunity-based rates
above.

Fully observed insertions must be in frame by default. The defaults are
research defaults, not universal or clinically validated thresholds.

## Normalization and consolidation

In repeated sequence, the same ALT can have multiple equivalent alignment-gap
placements. ITDiscover normalizes these to one canonical allele before fragment
counting. This is a representation correction, not consolidation.

With `--consolidate-minor-itd-variants`, an already-detected minor allele may
be assigned directly to one stronger allele only when their complete ALT
sequences, breakpoint shift, and support satisfy the configured consolidation
rules. A nearby, non-equivalent breakpoint can represent a distinct biological
allele or subclone, so consolidation is an optional error-suppression heuristic,
not proof of biological equivalence. Assignments are not chained, fragment
consensus is recalculated after assignment, and the TSV records every merge for
review.

See the [alignment walkthrough](alignment-walkthrough.md) for a concrete
fragment-by-fragment example.
