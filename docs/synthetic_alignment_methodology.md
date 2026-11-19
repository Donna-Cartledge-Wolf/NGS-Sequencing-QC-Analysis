# Version 1.2: ungapped synthetic nucleotide reference comparisons

## Purpose
Demonstrate how an apparently perfect base-level match can nevertheless be
ambiguous if the same read aligns equally well to two different references.
All references and reads are **invented random A/C/G/T strings** generated for
this portfolio exercise; labels *TOY_VIRAL_LIKE* and *TOY_HOST_LIKE* are purely
illustrative, **not** viral or human sequences.

## Algorithm
1. Parse Phred+33 FASTQ using the existing project's strict reader.
2. Parse FASTA references, rejecting duplicate labels and non-ACGT bases.
3. For each read, compare every full-length, ungapped placement on both the
   forward and reverse-complement strands of each toy reference.
4. Count mismatches and calculate percent identity as
   `100 * (read_length - mismatches) / read_length`.
5. Retain every equally best placement. Report a unique best placement only if
   exactly one placement is best; flag ties within a reference or across refs.
6. Apply illustrative maximum-mismatch and minimum-length rules.

## Interpretation
- A `CROSS_REFERENCE_AMBIGUOUS` read matches two different toy references equally
  well; **do not assign it as viral**.
- `UNIQUE_BEST_IN_TOY_REFERENCES` means unique only within this artificially tiny
  reference collection, **not** evidence of a biological virus.
- This project **does not calculate calibrated MAPQ**. `MAPQ=0` is often
  associated with uncertainty in common aligners, but no MAPQ should be
  fabricated from this exercise; output correctly says `NOT_CALCULATED`.
- Mean identity is not taxonomic confidence or assay sensitivity.

## Explicit non-features
No real virus/human sequence databases; no BLAST, Smith-Waterman/BWA/minimap2,
indels, gaps, soft clipping, scoring matrices, paired ends, base-quality-aware
alignment, adapter trimming, read duplicates, taxonomic classification, target
calling, biological inference, sample QC, or GMP/clinical validation.
The brute-force algorithm is deliberately not scalable and unsuitable for real
high-throughput NGS datasets.
