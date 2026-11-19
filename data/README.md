# Synthetic data only

`synthetic_reads.fastq` contains 12 artificially generated, 60-base nucleotide
reads originally prepared for an NGS quality-control teaching exercise.
Their sequence strings are synthetic and not derived from a biological sample,
viral experiment, patient, or client. Each FASTQ quality character is Phred+33.

Case studies include:

| Read | Teaching point |
|---|---|
| SIM_02 | 48 bases at Q35, then 12 bases at Q10; mean Q30 hides poor 3' tail. |
| SIM_03 | 60 bases at Q30; same mean Q30, much lower expected error burden. |
| SIM_05 | Poor-quality 5' segment; terminal 3' trimming alone does not fix it. |
| SIM_06 | Internal low-quality segment; terminal 3' trimming alone does not fix it. |
| SIM_09 | Ambiguous N calls; high mean quality does not eliminate unresolved bases. |

Do not upload real patient, sponsor, or confidential laboratory data to a public repository.


`synthetic_gapped_pairs.csv` contains only invented A/C/G/T strings for v1.3
Needleman-Wunsch and Smith-Waterman comparisons.
