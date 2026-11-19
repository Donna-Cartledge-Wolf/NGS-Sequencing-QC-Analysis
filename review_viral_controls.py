"""CLI for hypothetical, preclassified viral alignment-summary QC review."""
import argparse
import json
from ngs_qc.control_review import run

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='data/synthetic_control_summary.csv')
    parser.add_argument('--outdir', default='results')
    args = parser.parse_args()
    try:
        summary = run(args.input, args.outdir)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({k: v for k, v in summary.items() if k != 'results'}, indent=2))
    print(f"Review outputs written to: {args.outdir}")

if __name__ == '__main__':
    main()
