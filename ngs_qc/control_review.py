"""Educational, control-aware metrics for PRECLASSIFIED synthetic read counts.

No FASTQ alignment, taxonomic classification, reference search, viral
detection, limit-of-detection determination, or reportable call occurs.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

FIELDS = (
    'batch_id', 'sample_id', 'sample_type', 'virus_label',
    'total_filtered_reads', 'putative_viral_reads',
    'distinct_reference_regions', 'unique_supporting_reads',
    'mean_alignment_identity_pct'
)
ROLES = {'sample', 'extraction_blank', 'positive_process_control'}


@dataclass(frozen=True)
class ControlRecord:
    batch_id: str
    sample_id: str
    sample_type: str
    virus_label: str
    total_filtered_reads: int
    putative_viral_reads: int
    distinct_reference_regions: int
    unique_supporting_reads: int
    mean_alignment_identity_pct: float | None

    @property
    def rpm(self) -> float:
        return self.putative_viral_reads * 1_000_000 / self.total_filtered_reads


def load_records(path: str | Path) -> list[ControlRecord]:
    """Read hypothetical per-sample read-count summaries with basic validation."""
    with Path(path).open(newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or set(FIELDS) - set(reader.fieldnames):
            raise ValueError('CSV must include all documented columns')
        rows = list(reader)
    if not rows:
        raise ValueError('Input CSV has no records')
    records = []
    seen = set()
    depths = {}
    for idx, raw in enumerate(rows, start=2):
        def required(name: str) -> str:
            val = (raw[name] or '').strip()
            if not val:
                raise ValueError(f'Row {idx}: {name} cannot be empty')
            return val
        def number(name: str) -> int:
            val = required(name)
            try:
                parsed = int(val)
            except ValueError as exc:
                raise ValueError(f'Row {idx}: {name} must be an integer') from exc
            if parsed < 0:
                raise ValueError(f'Row {idx}: {name} cannot be negative')
            return parsed
        batch, sample, virus = (required('batch_id'), required('sample_id'), required('virus_label'))
        role = required('sample_type')
        if role not in ROLES:
            raise ValueError(f'Row {idx}: unrecognized sample_type {role!r}')
        key = (batch, sample, virus)
        if key in seen:
            raise ValueError(f'Row {idx}: duplicate batch/sample/virus combination')
        seen.add(key)
        depth = number('total_filtered_reads')
        hits = number('putative_viral_reads')
        regions = number('distinct_reference_regions')
        unique = number('unique_supporting_reads')
        if depth <= 0 or hits > depth or unique > hits or regions > hits:
            raise ValueError(f'Row {idx}: impossible depth/read/region counts')
        identity_text = (raw['mean_alignment_identity_pct'] or '').strip()
        try:
            identity = float(identity_text) if identity_text else None
        except ValueError as exc:
            raise ValueError(f'Row {idx}: malformed percent identity') from exc
        if hits and (identity is None or not (0 <= identity <= 100)):
            raise ValueError(f'Row {idx}: putative hits require a 0–100 identity')
        if not hits and identity is not None:
            raise ValueError(f'Row {idx}: identity requires at least one putative hit')
        sample_key = (batch, sample)
        if sample_key in depths and depths[sample_key] != (depth, role):
            raise ValueError(f'Row {idx}: inconsistent sample metadata across virus labels')
        depths[sample_key] = (depth, role)
        records.append(ControlRecord(batch, sample, role, virus, depth, hits,
                                     regions, unique, identity))
    return records


def make_review(records: list[ControlRecord]) -> list[dict]:
    """Compare sample vs same-batch, same-target extraction blank only.

    An undefined denominator stays undefined. No RPM ratio threshold is
    converted into positive/negative or any clinical decision.
    """
    blanks = {}
    for row in records:
        if row.sample_type == 'extraction_blank':
            key = (row.batch_id, row.virus_label)
            if key in blanks:
                raise ValueError(f'Duplicate extraction blanks for batch/target {key}')
            blanks[key] = row
    output = []
    for row in records:
        if row.sample_type != 'sample':
            continue
        blank = blanks.get((row.batch_id, row.virus_label))
        if blank is None:
            raise ValueError(f'No extraction blank for {row.batch_id}/{row.virus_label}')
        ratio = row.rpm / blank.rpm if blank.putative_viral_reads else None
        # Interpretation deliberately conservative, never a formal detection status.
        comments = []
        if row.putative_viral_reads:
            comments.append('Putative reads observed: inspect raw alignments and alternative references.')
        else:
            comments.append('No putative reads in supplied summary: does not demonstrate virus absence.')
        if blank.putative_viral_reads:
            comments.append('Extraction blank has putative reads: background/source requires review.')
        else:
            comments.append('Blank has zero summarized hits: ratio undefined, not an infinite signal.')
        if row.putative_viral_reads and row.unique_supporting_reads < row.putative_viral_reads:
            comments.append('Not all putative reads are uniquely supportive: examine multimapping.')
        if row.distinct_reference_regions > 1:
            comments.append('Multiple reported regions: verify independent, nonoverlapping mappings.')
        output.append({
            'batch_id': row.batch_id,
            'virus_label': row.virus_label,
            'sample_id': row.sample_id,
            'sample_depth': row.total_filtered_reads,
            'sample_putative_reads': row.putative_viral_reads,
            'sample_rpm': round(row.rpm, 6),
            'blank_id': blank.sample_id,
            'blank_depth': blank.total_filtered_reads,
            'blank_putative_reads': blank.putative_viral_reads,
            'blank_rpm': round(blank.rpm, 6),
            'sample_to_blank_rpm_ratio': round(ratio, 6) if ratio is not None else None,
            'distinct_reference_regions': row.distinct_reference_regions,
            'unique_supporting_reads': row.unique_supporting_reads,
            'mean_alignment_identity_pct': row.mean_alignment_identity_pct,
            'review_status': 'TECHNICAL_REVIEW_REQUIRED',
            'review_notes': ' '.join(comments),
        })
    if not output:
        raise ValueError('No client samples in input')
    return sorted(output, key=lambda x: (x['batch_id'], x['virus_label'], x['sample_id']))


def run(input_csv: str | Path, output_dir: str | Path) -> dict:
    rows = make_review(load_records(input_csv))
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / 'viral_control_review.csv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        'dataset': 'HYPOTHETICAL SYNTHETIC SUMMARY COUNTS - FOR EDUCATION ONLY',
        'sample_target_comparisons': len(rows),
        'extraction_blank_has_reads': sum(x['blank_putative_reads'] > 0 for x in rows),
        'ratio_undefined_blank_zero': sum(x['sample_to_blank_rpm_ratio'] is None for x in rows),
        'rule': 'All comparisons are technical review examples; NO positive/negative decisions',
        'limitation': ('Numbers represent PRECLASSIFIED fictional alignment summaries, not '
                       'independent FASTQ alignment, BLAST, pathogen detection, or a validated workflow.'),
        'results': rows,
    }
    (outdir / 'viral_control_review.json').write_text(
        json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary
