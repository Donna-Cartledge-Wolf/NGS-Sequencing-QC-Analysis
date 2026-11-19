# Development roadmap (post v1.3)

Completed educational modules:
- v1.0: synthetic FASTQ QC, Phred scoring and simple terminal trimming.
- v1.1: control-aware review of fictional, preclassified viral hit summaries.
- v1.2: ungapped full-read comparison of invented sequences against toy refs.
- v1.3: global Needleman-Wunsch and local Smith-Waterman gapped alignments of
  short fictional DNA pairs with linear gap penalties and basic scoring.

Potential next steps (not implemented):
1. Add affine gap costs and compare alternative traceback paths.
2. Adapter-aware trimming and representative quality profiles.
3. SAM/BAM interpretation with established bioinformatics software.
4. Database provenance, host subtraction, reference version control and
   decoy references within appropriate validated workflows.
5. Synthetic truth-set precision/recall, edge-case evaluation, and reproducible
   workflow packaging.

No real pathogen detection, taxonomic assignment, clinical reporting or GMP
validation is implied by these educational demonstrations.
