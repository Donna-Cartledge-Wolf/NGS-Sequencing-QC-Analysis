"""Transparent, UNGAPPED exact-length DNA read matching for a synthetic teaching dataset.

Not a production aligner. Does not infer taxonomy, estimate MAPQ, call variants,
or perform gapped/soft-clipped alignment. No biological reference sequences included.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import json
from .fastq import read_fastq

DNA = set("ACGT")
_RC = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def parse_fasta(path: str | Path) -> dict[str, str]:
    """Parse one-line or wrapped DNA FASTA, rejecting duplicate names and ambiguity."""
    references: dict[str, str] = {}
    name = None
    chunks = []
    def store():
        if name is not None:
            sequence = "".join(chunks).upper()
            if not sequence or set(sequence) - DNA:
                raise ValueError(f"Reference {name!r} must contain only A/C/G/T")
            references[name] = sequence
    with open(path, encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                name = line[1:].split()[0]
                if not name or name in references:
                    raise ValueError("Missing/duplicate reference identifier")
                chunks = []
            else:
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                chunks.append(line)
    store()
    if not references:
        raise ValueError("No references in FASTA")
    return references


@dataclass(frozen=True)
class Placement:
    reference: str
    position_1based: int
    strand: str
    mismatches: int


def all_placements(seq: str, references: dict[str, str]) -> list[Placement]:
    """Enumerate ungapped, full-length placements on both strands (no clipping)."""
    seq = seq.upper()
    if not seq or set(seq) - DNA:
        raise ValueError("This educational matcher requires unambiguous A/C/G/T reads")
    orientations = [("+", seq), ("-", reverse_complement(seq))]
    matches: list[Placement] = []
    for refname, reference in references.items():
        for strand, query in orientations:
            for i in range(len(reference) - len(query) + 1):
                mismatches = sum(a != b for a, b in zip(reference[i:i+len(query)], query))
                matches.append(Placement(refname, i + 1, strand, mismatches))
    return matches


def evaluate_read(name: str, seq: str, references: dict[str, str], max_mismatches: int = 2, min_length: int = 24) -> dict:
    if min_length < 1 or max_mismatches < 0:
        raise ValueError("Invalid thresholds")
    result = {"read_id": name, "read_length": len(seq), "best_mismatches": "",
              "best_identity_pct": "", "best_reference": "", "best_position_1based": "",
              "best_strand": "", "equally_best_placements": 0, "equally_best_references": "",
              "assignment": "", "review_note": "", "mapq": "NOT_CALCULATED"}
    if len(seq) < min_length:
        result.update(assignment="TOO_SHORT", review_note="Read below illustrative minimum length")
        return result
    if set(seq.upper()) - DNA:
        result.update(assignment="NON_ACGT", review_note="Ambiguous bases not handled in simplified matcher")
        return result
    placements = all_placements(seq, references)
    if not placements:
        result.update(assignment="NO_PLACEMENT", review_note="Read longer than every reference")
        return result
    least = min(item.mismatches for item in placements)
    if least > max_mismatches:
        result.update(assignment="NO_MATCH_WITHIN_LIMIT", best_mismatches=least,
                      review_note="No full-length alignment within illustrative mismatch limit")
        return result
    best = [item for item in placements if item.mismatches == least]
    refnames = sorted({hit.reference for hit in best})
    chosen = sorted(best, key=lambda hit: (hit.reference, hit.position_1based, hit.strand))[0]
    result.update(best_mismatches=least, best_identity_pct=round(100*(len(seq)-least)/len(seq), 3),
                  best_reference=chosen.reference if len(refnames)==1 else "AMBIGUOUS",
                  best_position_1based=chosen.position_1based if len(best)==1 else "",
                  best_strand=chosen.strand if len(best)==1 else "",
                  equally_best_placements=len(best), equally_best_references=";".join(refnames))
    if len(refnames)>1:
        result.update(assignment="CROSS_REFERENCE_AMBIGUOUS", review_note="Equal-best placements across references; do not infer source")
    elif len(best)>1:
        result.update(assignment="MULTIPLE_POSITIONS", review_note="Multiple equally best placements within one reference")
    else:
        result.update(assignment="UNIQUE_BEST_IN_TOY_REFERENCES", review_note="Unique only within tiny synthetic reference set; not evidence of biological origin")
    return result


def run(input_fastq: str | Path, fasta: str | Path, output_dir: str | Path, *, max_mismatches: int=2, min_length: int=24) -> dict:
    references = parse_fasta(fasta)
    reads = list(read_fastq(input_fastq))
    if not reads:
        raise ValueError("Empty FASTQ")
    records = [evaluate_read(read.name, read.sequence, references, max_mismatches, min_length) for read in reads]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "synthetic_reference_matches.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    summary = {"dataset":"ENTIRELY SYNTHETIC / EDUCATION ONLY", "reads":len(reads),
               "reference_count":len(references), "max_mismatches":max_mismatches,
               "minimum_read_length":min_length,
               "assignment_counts":{k:sum(r["assignment"]==k for r in records) for k in sorted({r["assignment"] for r in records})},
               "limitations":"Ungapped full-length matches only; toy references; no calibrated MAPQ, taxonomy, pathogen detection, or regulatory results."}
    (output_dir / "synthetic_reference_summary.json").write_text(json.dumps(summary,indent=2)+"\n", encoding="utf-8")
    return summary
