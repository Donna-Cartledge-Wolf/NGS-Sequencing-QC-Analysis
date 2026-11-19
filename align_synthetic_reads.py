"""Run the fully fictional v1.2 sequence-matching demonstration."""
import argparse
import json
from ngs_qc.reference_demo import run

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/synthetic_alignment_reads.fastq")
    parser.add_argument("--references", default="data/synthetic_references.fasta")
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--max-mismatches", type=int, default=2)
    parser.add_argument("--min-length", type=int, default=24)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.references, args.outdir,
                         max_mismatches=args.max_mismatches,
                         min_length=args.min_length), indent=2))
