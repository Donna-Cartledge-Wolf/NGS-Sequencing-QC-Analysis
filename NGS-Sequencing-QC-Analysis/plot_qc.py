"""Optional plots: pip install -r requirements-visuals.txt; python plot_qc.py."""
import csv
from pathlib import Path


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = Path(__file__).resolve().parent
    results = root / "results"
    images = root / "assets"
    images.mkdir(exist_ok=True)
    with (results / "qc_by_position.csv").open(newline="", encoding="utf-8") as stream:
        positions = list(csv.DictReader(stream))
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    ax.plot([int(r["position_1based"]) for r in positions],
            [float(r["mean_phred_arithmetic"]) for r in positions],
            linewidth=2)
    ax.set(title="Mean base-call quality by read position (synthetic data)",
           xlabel="Position along read (1-based)",
           ylabel="Arithmetic mean Phred quality")
    ax.axhline(30, linestyle="--", linewidth=1, label="Q30 guide (illustrative)")
    ax.legend(loc="lower left")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(images / "quality_by_position.png", dpi=175)
    plt.close(fig)

    with (results / "qc_per_read.csv").open(newline="", encoding="utf-8") as stream:
        reads = list(csv.DictReader(stream))
    ids = [r["read_id"] for r in reads]
    pre = [float(r["expected_base_errors"]) for r in reads]
    post = [float(r["post_trim_expected_errors"]) if r["retained_after_min_length"] == "True" else 0
            for r in reads]
    x = list(range(len(ids)))
    fig, ax = plt.subplots(figsize=(10, 4.5))
    w = 0.36
    ax.bar([i - w/2 for i in x], pre, w, label="Before trimming")
    ax.bar([i + w/2 for i in x], post, w, label="After trimming")
    ax.set(title="Expected base-call errors per read (synthetic data)",
           ylabel="Sum of per-base error probabilities")
    ax.set_xticks(x, ids, rotation=45, ha="right")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(images / "expected_errors_before_after.png", dpi=175)
    plt.close(fig)
    print("Created assets/quality_by_position.png and assets/expected_errors_before_after.png")


if __name__ == "__main__":
    main()
