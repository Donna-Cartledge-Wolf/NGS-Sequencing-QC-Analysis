"""Transparent educational read-level and position-level Phred statistics."""
import csv
import json
from collections import defaultdict
from math import prod
from pathlib import Path
from .fastq import Read, read_fastq

# Illustrative REVIEW triggers, deliberately not validated pass/fail criteria.
REVIEW_EXPECTED_ERRORS = 0.25
REVIEW_Q30_PERCENT = 80.0
REVIEW_TAIL10_MEAN_Q = 20.0


def error_probabilities(scores: list[int]) -> list[float]:
    """Calibrated Phred estimate: p(error) = 10**(-Q/10), per base."""
    return [10 ** (-q / 10.0) for q in scores]


def stats_for_read(read: Read) -> dict:
    qs = read.scores
    ps = error_probabilities(qs)
    q30 = sum(q >= 30 for q in qs)
    tail_q = qs[-10:]
    expected = sum(ps)
    flags = []
    if "N" in read.sequence:
        flags.append("contains_N")
    if 100 * q30 / len(qs) < REVIEW_Q30_PERCENT:
        flags.append("low_Q30_fraction")
    if sum(tail_q) / len(tail_q) < REVIEW_TAIL10_MEAN_Q:
        flags.append("low_3prime_tail_quality")
    if expected > REVIEW_EXPECTED_ERRORS:
        flags.append("elevated_expected_errors")
    return {
        "read_id": read.name,
        "read_length": len(qs),
        "mean_phred_arithmetic": round(sum(qs) / len(qs), 3),
        "pct_bases_Q30_or_above": round(100 * q30 / len(qs), 3),
        "tail_last10_mean_phred": round(sum(tail_q) / len(tail_q), 3),
        "expected_base_errors": round(expected, 6),
        "probability_ge1_error_if_independent": round(1 - prod((1 - p) for p in ps), 6),
        "N_count": read.sequence.count("N"),
        "review_flags": ";".join(flags) if flags else "none",
    }


def trim_3prime_below_q(read: Read, q_threshold: int = 20) -> Read:
    """Remove only the contiguous 3' suffix with per-base quality below Q threshold.

    Does not trim adapters or internal low-quality positions. This simple rule
    is for illustrating the sensitivity / specificity tradeoff, not a choice of
    optimal settings for any real NGS platform or regulated assay.
    """
    scores = read.scores
    boundary = len(scores)
    while boundary > 0 and scores[boundary - 1] < q_threshold:
        boundary -= 1
    return Read(read.name, read.sequence[:boundary], read.encoded_quality[:boundary])


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run(input_fastq: str | Path, output_dir: str | Path, trim_q: int = 20,
        min_length: int = 35) -> dict:
    """Produce deterministic CSV, JSON, and retained trimmed FASTQ outputs."""
    if not 0 <= trim_q <= 93:
        raise ValueError("--trim-q must be between 0 and 93")
    if min_length < 1:
        raise ValueError("--min-length must be at least 1")
    reads = list(read_fastq(input_fastq))
    if not reads:
        raise ValueError("The FASTQ file contains no reads")
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)

    total_bases = 0
    raw_q30 = 0
    per_read = []
    by_position = defaultdict(list)
    retained = []
    trimmed_bases = 0
    total_expected_pre = 0.0
    total_expected_post = 0.0

    for read in reads:
        metrics = stats_for_read(read)
        qs = read.scores
        total_bases += len(qs)
        raw_q30 += sum(q >= 30 for q in qs)
        total_expected_pre += sum(error_probabilities(qs))
        for position, q in enumerate(qs, start=1):
            by_position[position].append(q)

        short = trim_3prime_below_q(read, q_threshold=trim_q)
        trim_count = len(read.sequence) - len(short.sequence)
        trimmed_bases += trim_count
        keep = len(short.sequence) >= min_length
        metrics.update({
            "trimmed_bases_3prime": trim_count,
            "post_trim_length": len(short.sequence),
            "post_trim_expected_errors": round(
                sum(error_probabilities(short.scores)), 6),
            "retained_after_min_length": keep,
        })
        per_read.append(metrics)
        if keep:
            retained.append(short)
            total_expected_post += sum(error_probabilities(short.scores))

    by_position_rows = [{
        "position_1based": pos,
        "num_reads_covering_position": len(qs),
        "mean_phred_arithmetic": round(sum(qs) / len(qs), 3),
        "pct_bases_Q30_or_above": round(
            100 * sum(q >= 30 for q in qs) / len(qs), 3),
    } for pos, qs in sorted(by_position.items())]

    summary = {
        "dataset": "SYNTHETIC / FOR EDUCATION ONLY",
        "input_file": Path(input_fastq).name,
        "input_reads": len(reads),
        "input_bases": total_bases,
        "Q30_or_above_bases": raw_q30,
        "pct_Q30_or_above_bases": round(100 * raw_q30 / total_bases, 3),
        "reads_flagged_for_review": sum(
            1 for row in per_read if row["review_flags"] != "none"),
        "example_review_rule": {
            "min_Q30_percentage": REVIEW_Q30_PERCENT,
            "max_expected_errors": REVIEW_EXPECTED_ERRORS,
            "min_tail10_mean_Q": REVIEW_TAIL10_MEAN_Q,
            "note": "Instructional flags only, not a validated pass/fail decision",
        },
        "trim_3prime_below_phred": trim_q,
        "minimum_retained_length_bases": min_length,
        "reads_retained_after_trimming": len(retained),
        "reads_excluded_for_minimum_length": len(reads) - len(retained),
        "bases_removed_from_3prime_tails": trimmed_bases,
        "total_expected_errors_before_trimming": round(total_expected_pre, 6),
        "total_expected_errors_across_retained_reads_after_trimming": round(
            total_expected_post, 6),
        "interpretation_warning": (
            "No read alignment, viral classification, contamination assessment, "
            "diagnostic result, or regulatory validation is performed."
        ),
    }
    write_csv(folder / "qc_per_read.csv", per_read)
    write_csv(folder / "qc_by_position.csv", by_position_rows)
    (folder / "qc_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (folder / "trimmed_reads.fastq").write_text(
        "".join(read.to_fastq() for read in retained), encoding="ascii")
    return summary
