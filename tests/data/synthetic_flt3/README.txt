Reference: reference.fa

This fixture set is designed to produce calls with ITDiscover's default
parameters.

The synthetic reads are built from the bundled 329-base FLT3 amplicon
reference sequence, so the following primer sequences can be used to test
trimming:
  forward-primer = GCAATTTAGGTATGAAAGCCAGC
  reverse-primer = CTTTCAGCATTTTGACGGCAACC

Fragments:
  wt-1..4             = wild type controls
  exact-itd-1,4,5     = exact 15 bp ITD copied from zero-based reference
                        positions 79-93.
  fuzzy-itd-1         = 6 bp insertion with one mismatch in the copied segment;
                        it should appear with
                        `--max-copy-mismatch-rate 0.166667`.
  fuzzy-itd-2         = 6 bp insertion with two mismatches in the copied
                        segment; it should only appear with
                        `--max-copy-mismatch-rate 0.333334`.
  spacer-itd-1..2     = insertion with a 10 bp copied segment from reference
                        positions 79-88 and spacers of three bases on the left
                        and two on the right (AAA + copied segment + TT).
  trailing-ins-1      = trailing insertion example; the inserted C-rich tail
                        is meant to illustrate a trailing insertion and is
                        not expected to call as an ITD.

With the primers above, default exact-mode calling reports:
  - one passing 15 bp ITD with 3 mutant and 9 wild-type fragments;
  - one spacer-containing candidate with 2 mutant and 10 wild-type fragments,
    filtered by LOW_MUTANT_FRAGMENT_COUNT.

The sample outcome is ITD detected. Its QC status is warn because the filtered
candidate is retained in the report.
