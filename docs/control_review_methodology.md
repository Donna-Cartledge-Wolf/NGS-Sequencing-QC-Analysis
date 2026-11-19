# Control-aware viral hit triage: Synthetic teaching demonstration

**Version 1.1 extension of the Wolf Analytics NGS QC project.**

This component compares **hypothetical, already classified read-count summaries**
from a client sample against the same-run, same-target extraction blank.
It does **not** analyze FASTQ sequences, perform sequence alignment, call
contaminants, or establish the absence of viruses.

## Input fields

- `batch_id`: simulated processing group.
- `sample_id` and `sample_type`: client sample, extraction blank, or positive process control.
- `virus_label`: fictional label, not a biological sequence.
- `total_filtered_reads`: number of passing reads supplied by an imagined upstream process.
- `putative_viral_reads`: fictional reads provisionally assigned to a virus.
- `distinct_reference_regions`: fictional count of potentially independent reference areas.
- `unique_supporting_reads`: fictional subset considered uniquely mapping by the imagined upstream step.
- `mean_alignment_identity_pct`: fictional upstream annotation; not sufficient alone for identification.

The program verifies basic numeric consistency and that an extraction blank
exists for each batch/virus label being reviewed.

## Core calculations

```text
RPM = putative_reads / total_filtered_reads * 1,000,000
sample-to-blank ratio = sample_RPM / blank_RPM     (only if blank_RPM > 0)
```

For **CLIENT_A, Virus_X**:

```text
Sample RPM = 7 / 1,200,000 * 1,000,000 = 5.833333
Blank RPM = 2 / 500,000 * 1,000,000 = 4.000000
Ratio = 5.833333 / 4.000000 = 1.458333
```

A 1.46 ratio does not define pass, fail, positive, or negative. With only a
few putative reads the ratio also has considerable sampling uncertainty.
Results must be evaluated against *assay-specific validated criteria*,
full controls, read uniqueness, genome breadth, database provenance, and
appropriately designed confirmation methods.

For **CLIENT_C, Virus_Y**, the extraction blank has zero putative reads.
Dividing by blank RPM is mathematically undefined: the script stores a
JSON `null` and an empty CSV field, **not** `Infinity`, a pseudo-count, or
an automatic positive result. Whether to use pseudocounts would be an
independent, explicitly justified analytical design decision.

## Scientific and quality caveats

A nonzero signal in the blank may indicate background contamination,
barcode misassignment, or annotation artifacts, but source attribution
cannot be made from a summary table. Even zero summary reads in a blank
cannot establish an absence of contamination or guarantee sensitivity.
Multiple reported genomic regions may offer stronger supporting evidence
**only after** inspecting the independent alignments, host homology,
sequence complexity, and database quality.

**No numerical RPM ratio or identity cutoff is universally applicable** to
all viral safety assays. No regulated QC release determination is made.
This software is educational and not GMP validated or for clinical use.

## Reproduction

From the repository root:

```cmd
py review_viral_controls.py --input data/synthetic_control_summary.csv --outdir results
py -m unittest discover -s tests -v
```

Output: `results/viral_control_review.csv` and
`results/viral_control_review.json`.

## Next possible development

A later, separate module could implement transparent *toy alignment* against
synthetic host-like/virus-like sequences, permitting inspection of ambiguous
mappings and mapping quality. That functionality is **not implemented here**.
