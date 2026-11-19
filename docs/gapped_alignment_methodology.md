# Version 1.3: Global and local gapped alignment of fictional DNA

## Purpose and input provenance

This module implements a small teaching example using invented DNA strings
in `data/synthetic_gapped_pairs.csv`. These strings are **not** actual human,
viral, or other biological sequences. Comparisons are not evidence of viral
presence, pathogenicity, a clinical result, or even reference uniqueness.

## Methods

- **Needleman-Wunsch (global):** dynamic programming aligns both complete
  sequences, including their terminal bases.
- **Smith-Waterman (local):** dynamic programming finds one highest-scoring
  matching region within the two sequences; unaligned flanks are omitted.
- Default scores are +2 for a match, -1 for a substitution, and -2 for each
  gap base. The gap penalty is **linear**, not affine (no separate gap-opening
  and gap-extension costs). Score choice can change the resulting alignment.
- A deterministic traceback prioritizes diagonal, deletion, then insertion
  when multiple moves have the same score. It does **not** enumerate all
  equally optimal alignments, a major limitation for repetitive sequences.
- Local matches tie by the first maximum visited in row-major order. A local
  score of zero is explicitly reported as `NO_POSITIVE_LOCAL_ALIGNMENT`.
- Only unambiguous A/C/G/T sequences of at most 500 bases are accepted.
  There is no read-quality weighting, reverse complement search, indexing,
  paired-read modeling, indel realignment around biological variants, or
  calibrated mapping quality (MAPQ).

## Interpretation of output

Positions are 1-based inclusive and refer to the *chosen* aligned region.
A local alignment may have unaligned terminal bases in both sequences.

Definitions relative to the **query** sequence (the second input):
- **Substitution (X):** two aligned nucleotides disagree.
- **Insertion (I):** a nucleotide appears in the query opposite a reference gap.
- **Deletion (D):** a reference nucleotide appears opposite a query gap.
- **Match (=):** identical aligned nucleotides.

The `operations` field resembles extended CIGAR with `=`, `X`, `I`, and `D`,
with repeated operations encoded by a count. It is a teaching annotation,
**not** a complete SAM/CIGAR file and has no soft-clipping operations.

**Percent identity denominator:** matching positions divided by the full
number of aligned columns **including gap columns**, times 100. This is one
explicit identity convention; some alignment software uses other denominators.
A 100% identity for a local alignment says nothing about its unaligned flanks.

## What the output does NOT establish

One optimal traceback does not prove unique gap placement. For example,
`HOMOPOLYMER_GAP` has many equivalent single-gap placements. These toy
alignments are not taxonomic classifications, true variant calls, clinical
NGS pipelines, regulated tests, or proof of a software validation state.
No real sample, proprietary data, or patient data should be uploaded.

## Run and independent checks

```cmd
py align_gapped_synthetic.py
py -m unittest discover -s tests -v
```

Read `results/synthetic_gapped_alignments.csv` in Excel. For the substitution
example, verify 7 of 8 bases match (87.5% identity). For the insertion example,
verify 8 matching columns and 1 insertion column (8/9 * 100 = 88.8889%).
For the deletion example, verify 7 matching columns and 1 deletion (7/8 * 100
= 87.5%). Re-check the chosen alignment strings before interpreting numbers.

**Scoring sensitivity experiment:** rerun with a different gap penalty and
compare output in a separate folder, so original example results remain intact:

```cmd
py align_gapped_synthetic.py --gap -5 --outdir results/scoring_gap_minus5
```

For `SCORING_SENSITIVITY` (`AAA` versus `AAC`), the default scheme favors
`2=1X`: two matching bases and one mismatch, total score 3. With `--mismatch -3
--gap -1`, an alignment using one query insertion and one query deletion can
score 2, beating the two-match/one-mismatch score of 1. This illustrates that
reported gap placements and event categories can change solely with scoring.

Results under alternative scoring do not imply biological correctness.
