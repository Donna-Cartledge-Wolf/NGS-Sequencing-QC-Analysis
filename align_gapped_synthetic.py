"""Run teaching-only global/local alignments of wholly fictional DNA pairs."""
import argparse
import json
from ngs_qc.gapped_alignment import Scoring, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/synthetic_gapped_pairs.csv")
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--match", type=int, default=2)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-2)
    args = parser.parse_args()
    scoring = Scoring(args.match, args.mismatch, args.gap)
    print(json.dumps(run(args.input, args.outdir, scoring), indent=2))


if __name__ == "__main__":
    main()
