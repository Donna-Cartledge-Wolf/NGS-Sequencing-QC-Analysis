# Planned follow-on development (beyond v1.1)

1. **Adapter identification and trimming** with platform-appropriate logic.
2. **Optional alignment exercise** against synthetic human-like and viral-like
   reference sequences; demonstrate ambiguous multi-mapping and MAPQ 0.
3. **Extend synthetic control review (v1.1 added)**: the current module compares
   hypothetical normalized counts, blanks, and supplied region summaries.
   Future work: examine actual synthetic read alignments and evidence for or
   against assignments.
4. **Reference database and pipeline traceability** (version pins, provenance,
   reference indexes, test cases, checksums, and per-run metadata).
5. **More complete validation-style assessment** using synthetic truth sets and
   sensitivity/specificity measures, documented failure modes, and reproducible
   parameter selection.

Nothing in this roadmap should be described as existing functionality, a
clinical application, or a regulated NGS test.
