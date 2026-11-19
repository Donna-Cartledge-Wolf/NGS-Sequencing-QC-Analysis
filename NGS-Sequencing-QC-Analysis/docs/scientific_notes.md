# Methodology and interpretation notes

## Phred mathematics

For a base-call quality score **Q**:

- Estimated probability of an incorrect base: **p = 10^(-Q/10)**.
- **Expected number of errors per read**: `sum(p_i)` across base positions.
- **Probability of >=1 error**, assuming *independent* and well-calibrated errors:
  `1 - product(1 - p_i)`.

These are mathematical estimates; Phred scores do not independently validate
biological provenance, correct sample identification, or error-free library
preparation. Real base errors may be correlated, and quality calibration can vary.

## Why arithmetic mean Q alone can mislead

Both SIM_02 and SIM_03 have arithmetic mean Q30:

- SIM_02: `48 x Q35 + 12 x Q10` -> mean Q30 and **1.215179** expected errors.
- SIM_03: `60 x Q30` -> mean Q30 and **0.060000** expected errors.

Mean Q is the arithmetic mean of logarithmic scores, **not** the error
probability of an 'average' base. Interpret Q distribution and expected errors.

## Trimming method in v1

The demonstration removes contiguous low-quality bases from the **3' end only**
until the terminal remaining base has Q >= the chosen threshold (default Q20).
A read is retained only if it meets the illustrative minimum length (35 nt).
This intentionally **does not** remove adapter sequences, internal low-quality
bases, or low-quality 5' sequence. It is not equivalent to optimized `fastp`,
`Trimmomatic`, or alignment-aware filtering. The script preserves the original
sample-quality symbols in the trimmed FASTQ. Review flags are educational;
there are no diagnostic or GMP accept/reject criteria here.

## Limitations

- Requires standard four-line FASTQ records with Phred+33 ASCII encoding;
  wrapped/multiline records and Phred+64 are unsupported.
- Input permits IUPAC ambiguous nucleotide symbols; only literal `N` is
  counted as ambiguous in the summary metric.
- Independent-error probability is only a model; use caution for correlated
  errors such as homopolymer artifacts.
- Reads are **not** aligned, assigned taxonomy, or assessed for viruses.
- Quality metrics alone **cannot** establish the presence/absence of viruses.
- The synthetic input does not assess sensitivity, specificity, detection
  limits, reproducibility, or performance on real samples.

## Relevance to laboratory quality systems

This project demonstrates quality-metric calculations, transparent parameter
selection, example-data provenance, automated tests, and reproducible outputs.
It **does not** demonstrate regulatory qualification/validation of the
underlying NGS assay or bioinformatics pipeline. Actual GMP workflows require
approved SOPs, validation to the intended use, data-integrity controls,
traceability, change control, and appropriate quality oversight.
