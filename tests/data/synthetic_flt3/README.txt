Reference: reference.fasta

This fixture set is designed to produce calls with ITDiscover's default
parameters.

The reads are built from the full FLT3 reference so the following primer
sequences can be used to test trimming:
  forward-primer = GCAATTTAGGTATGAAAGCCAGCTAC
  reverse-primer = CTTTCAGCATTTTGACGGCAACC

Fragments:
  wt-1..4             = wild type controls
  exact-itd-1..5      = exact 15 bp ITD copied from reference positions
                        80-94.
  fuzzy-itd-1         = 6 bp insertion with one mismatch in the copied tract;
                        it should appear with `--max-mismatches 1`.
  fuzzy-itd-2         = 6 bp insertion with two mismatches in the copied
                        tract; it should only appear with `--max-mismatches 2`.
  spacer-itd-1..2     = insertion with a 10 bp copied tract from reference
                        positions 80-89 and spacers of three bases on the
                        left and two on the right (AAA + copied tract + TT).
  trailing-ins-1      = trailing insertion example; the inserted C-rich tail
                        is meant to illustrate a trailing insertion and is
                        not expected to call as an ITD.

Default exact-mode calling should report two ITDs: the 15 bp exact duplication
and the spacer-containing duplication.
