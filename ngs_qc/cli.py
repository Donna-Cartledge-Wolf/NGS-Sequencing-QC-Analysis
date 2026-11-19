"""Command line entry point."""
import argparse
import json
from pathlib import Path
from .analysis import run


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Educational FASTQ quality-control demo using synthetic data.")
    parser.add_argument("--input", type=Path,
                        default=root / "data" / "synthetic_reads.fastq")
    parser.add_argument("--outdir", type=Path, default=root / "results")
    parser.add_argument("--trim-q", type=int, default=20,
                        help="trim contiguous terminal 3' bases below this Q (default: 20)")
    parser.add_argument("--min-length", type=int, default=35,
                        help="minimum length to retain a trimmed read (default: 35)")
    args = parser.parse_args()
    try:
        report = run(args.input, args.outdir, args.trim_q, args.min_length)
    except (ValueError, OSError, UnicodeError) as exc:
        parser.exit(status=2, message=f"FASTQ QC error: {exc}\n")
    print(json.dumps(report, indent=2))
    print(f"\nResults written to: {args.outdir.resolve()}")
