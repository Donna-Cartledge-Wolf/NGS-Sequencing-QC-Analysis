"""Small linear-gap global/local aligner for non-biological teaching sequences.

Educational algorithms: Needleman-Wunsch (global) and Smith-Waterman (local).
Not appropriate for NGS mapping, taxonomy, mutation reporting or clinical use.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import json

DNA = frozenset("ACGT")
MAX_SEQUENCE_LENGTH = 500  # Bound time and memory for this teaching-only DP.


@dataclass(frozen=True)
class Scoring:
    match: int = 2
    mismatch: int = -1
    gap: int = -2

    def __post_init__(self):
        if not all(type(x) is int for x in (self.match, self.mismatch, self.gap)):
            raise ValueError("Scoring parameters must be integers")
        if self.match <= 0 or self.mismatch > 0 or self.gap >= 0:
            raise ValueError("Require match > 0, mismatch <= 0, gap < 0")


def _validate_sequence(seq: str, field: str) -> str:
    if not isinstance(seq, str):
        raise ValueError(f"{field} must be a string")
    seq = seq.upper()
    if not seq or len(seq) > MAX_SEQUENCE_LENGTH or set(seq) - DNA:
        raise ValueError(f"{field} must be 1 to {MAX_SEQUENCE_LENGTH} A/C/G/T bases")
    return seq


def _rle(ops: list[str]) -> str:
    """Extended CIGAR-like notation: = identical, X substitution, I, D."""
    if not ops:
        return ""
    result = []
    last = ops[0]
    count = 0
    for op in ops:
        if op != last:
            result.append(f"{count}{last}")
            last = op
            count = 0
        count += 1
    result.append(f"{count}{last}")
    return "".join(result)


@dataclass(frozen=True)
class Alignment:
    method: str
    score: int
    aligned_reference: str
    aligned_query: str
    reference_start_1based: int | None
    reference_end_1based: int | None
    query_start_1based: int | None
    query_end_1based: int | None
    matches: int
    substitutions: int
    insertions_in_query: int
    deletions_in_query: int
    identity_pct: float | None
    operations: str
    status: str

    def as_dict(self) -> dict:
        return vars(self).copy()


def align(reference: str, query: str, method: str = "global", scoring: Scoring | None = None) -> Alignment:
    """Align complete synthetic sequences, returning a deterministic optimum.

    Scores use a linear gap penalty per deleted/inserted base, not affine gaps.
    Tied optimal paths are not enumerated, and no mapping uniqueness is implied.
    Rows consume reference bases; columns consume query bases.
    """
    reference = _validate_sequence(reference, "reference")
    query = _validate_sequence(query, "query")
    if method not in {"global", "local"}:
        raise ValueError("method must be 'global' or 'local'")
    scoring = scoring if scoring is not None else Scoring()
    if not isinstance(scoring, Scoring):
        raise ValueError("scoring must be a Scoring instance")
    n, m = len(reference), len(query)
    dp = [[0] * (m+1) for _ in range(n+1)]
    trace = [[""] * (m+1) for _ in range(n+1)]
    if method == "global":
        for i in range(1,n+1):
            dp[i][0] = i * scoring.gap
            trace[i][0] = "D"  # gap in query
        for j in range(1,m+1):
            dp[0][j] = j * scoring.gap
            trace[0][j] = "I"  # gap in reference
    best_score, best_i, best_j = 0, 0, 0
    for i in range(1,n+1):
        for j in range(1,m+1):
            diagonal = dp[i-1][j-1] + (scoring.match if reference[i-1] == query[j-1] else scoring.mismatch)
            deletion = dp[i-1][j] + scoring.gap
            insertion = dp[i][j-1] + scoring.gap
            # Deterministic tie priority: diagonal, deletion, insertion.
            candidates = [(diagonal, "M"), (deletion, "D"), (insertion, "I")]
            optimum = max(x[0] for x in candidates)
            if method == "local" and optimum <= 0:
                dp[i][j], trace[i][j] = 0, ""
            else:
                dp[i][j] = optimum
                trace[i][j] = next(op for score,op in candidates if score == optimum)
            # Choose the first row-major local maximum if multiple equal maxima.
            if method == "local" and dp[i][j] > best_score:
                best_score, best_i, best_j = dp[i][j], i, j
    if method == "global":
        i, j, score = n, m, dp[n][m]
    else:
        i, j, score = best_i, best_j, best_score
    end_i, end_j = i, j
    if method == "local" and score == 0:
        return Alignment(method, 0, "", "", None, None, None, None,
                         0, 0, 0, 0, None, "", "NO_POSITIVE_LOCAL_ALIGNMENT")

    aligned_ref_rev, aligned_qry_rev, operations_rev = [], [], []
    while (i > 0 or j > 0) if method == "global" else dp[i][j] > 0:
        action = trace[i][j]
        if action == "M":
            rbase, qbase = reference[i-1], query[j-1]
            aligned_ref_rev.append(rbase)
            aligned_qry_rev.append(qbase)
            operations_rev.append("=" if rbase == qbase else "X")
            i -= 1
            j -= 1
        elif action == "D":
            aligned_ref_rev.append(reference[i-1])
            aligned_qry_rev.append("-")
            operations_rev.append("D")
            i -= 1
        elif action == "I":
            aligned_ref_rev.append("-")
            aligned_qry_rev.append(query[j-1])
            operations_rev.append("I")
            j -= 1
        else:
            raise RuntimeError("Invalid traceback state")
    operations = list(reversed(operations_rev))
    matched = operations.count("=")
    substitution = operations.count("X")
    insertion = operations.count("I")
    deletion = operations.count("D")
    aligned_length = len(operations)
    return Alignment(
        method=method, score=score,
        aligned_reference="".join(reversed(aligned_ref_rev)),
        aligned_query="".join(reversed(aligned_qry_rev)),
        reference_start_1based=i+1, reference_end_1based=end_i,
        query_start_1based=j+1, query_end_1based=end_j,
        matches=matched, substitutions=substitution,
        insertions_in_query=insertion, deletions_in_query=deletion,
        identity_pct=round(100*matched/aligned_length, 4),
        operations=_rle(operations), status="ALIGNED",
    )


def read_pairs(path: str | Path) -> list[dict]:
    """Read self-contained synthetic sequence pairs, rejecting duplicate IDs."""
    records = []
    used = set()
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"pair_id", "reference", "query", "teaching_focus"}
        if not reader.fieldnames or set(reader.fieldnames) != required:
            raise ValueError("Pair CSV must have pair_id, reference, query, teaching_focus columns")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed CSV row")
            name = row["pair_id"].strip()
            if not name or name in used:
                raise ValueError("Pair identifiers must be unique and non-empty")
            used.add(name)
            records.append({
                "pair_id": name,
                "reference": _validate_sequence(row["reference"].strip(), "reference"),
                "query": _validate_sequence(row["query"].strip(), "query"),
                "teaching_focus": row["teaching_focus"].strip(),
            })
    if not records:
        raise ValueError("Pair CSV contains no examples")
    return records


def run(input_csv: str | Path, outdir: str | Path, scoring: Scoring | None = None) -> dict:
    scoring = scoring or Scoring()
    pairs = read_pairs(input_csv)
    records = []
    for pair in pairs:
        for method in ("global", "local"):
            result = align(pair["reference"], pair["query"], method, scoring)
            records.append({"pair_id":pair["pair_id"], "teaching_focus":pair["teaching_focus"],
                            "reference":pair["reference"], "query":pair["query"],
                            **result.as_dict()})
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "synthetic_gapped_alignments.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    summary = {
        "dataset":"ENTIRELY SYNTHETIC / EDUCATION ONLY",
        "pair_count":len(pairs), "alignment_count":len(records),
        "scoring":vars(scoring),
        "methods":["Needleman-Wunsch (global)", "Smith-Waterman (local)"],
        "identity_denominator":"All alignment columns, including gaps",
        "cigar":"Extended teaching notation: = exact, X substitution, I query insertion, D query deletion",
        "limits":"Linear gap penalty; one deterministically chosen optimum; no affine gaps, uniqueness, MAPQ, biological sequence interpretation or clinical validation",
    }
    (output / "synthetic_gapped_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary
