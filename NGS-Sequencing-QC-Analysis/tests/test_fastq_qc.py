import gzip
import json
from pathlib import Path
import tempfile
import unittest

from ngs_qc.fastq import Read, read_fastq
from ngs_qc.analysis import (
    error_probabilities, stats_for_read, trim_3prime_below_q, run,
)

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "data" / "synthetic_reads.fastq"


class FastqTest(unittest.TestCase):
    def test_dataset_has_12_reads(self):
        self.assertEqual(len(list(read_fastq(EXAMPLE))), 12)

    def test_sim02_matches_training_exercise(self):
        read = list(read_fastq(EXAMPLE))[1]
        self.assertEqual(read.name, "SIM_02")
        self.assertEqual(read.scores, [35] * 48 + [10] * 12)
        self.assertEqual(stats_for_read(read)["mean_phred_arithmetic"], 30.0)
        self.assertAlmostEqual(sum(error_probabilities(read.scores)), 1.2151789327688083)

    def test_sim03_uniform_q30(self):
        read = list(read_fastq(EXAMPLE))[2]
        self.assertEqual(read.scores, [30] * 60)
        self.assertAlmostEqual(stats_for_read(read)["expected_base_errors"], 0.06)

    def test_trimming_removes_only_3prime_suffix(self):
        read = list(read_fastq(EXAMPLE))[1]
        trimmed = trim_3prime_below_q(read, q_threshold=20)
        self.assertEqual(len(trimmed.sequence), 48)
        self.assertEqual(trimmed.scores, [35] * 48)
        self.assertAlmostEqual(sum(error_probabilities(trimmed.scores)),
                               48 * 10 ** -3.5)

    def test_zero_trimming_when_tail_good(self):
        read = list(read_fastq(EXAMPLE))[2]
        self.assertEqual(trim_3prime_below_q(read), read)

    def test_internal_low_quality_not_fixed_by_tail_trimming(self):
        read = list(read_fastq(EXAMPLE))[5]  # SIM_06 has internal low-quality bases
        self.assertEqual(trim_3prime_below_q(read), read)

    def test_bad_fastq_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "broken.fastq"
            path.write_text("@bad\nACGT\n+\n!!!\n")
            with self.assertRaises(ValueError):
                list(read_fastq(path))

    def test_truncated_fastq_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "broken.fastq"
            path.write_text("@bad\nACGT\n+\n")
            with self.assertRaises(ValueError):
                list(read_fastq(path))

    def test_gzip_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq.gz"
            with gzip.open(path, "wt", encoding="ascii") as stream:
                stream.write("@one\nACGT\n+\nIIII\n")
            self.assertEqual(len(list(read_fastq(path))), 1)

    def test_empty_fastq_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.fastq"
            path.write_text("")
            with self.assertRaisesRegex(ValueError, "no reads"):
                run(path, Path(tmp) / "out")

    def test_end_to_end_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "output"
            summary = run(EXAMPLE, out)
            self.assertEqual(summary["input_reads"], 12)
            self.assertEqual(summary["input_bases"], 720)
            self.assertEqual(summary["Q30_or_above_bases"], 611)
            self.assertEqual(summary["reads_flagged_for_review"], 6)
            self.assertEqual(summary["bases_removed_from_3prime_tails"], 14)
            self.assertEqual(summary["reads_retained_after_trimming"], 12)
            self.assertTrue((out / "trimmed_reads.fastq").is_file())
            self.assertEqual(len(list(read_fastq(out / "trimmed_reads.fastq"))), 12)
            self.assertEqual(json.loads((out / "qc_summary.json").read_text())["input_reads"], 12)

    def test_min_length_excludes_short_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = run(EXAMPLE, Path(tmp) / "out", min_length=55)
            self.assertEqual(summary["reads_retained_after_trimming"], 11)

    def test_reject_invalid_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                run(EXAMPLE, tmp, trim_q=-1)
            with self.assertRaises(ValueError):
                run(EXAMPLE, tmp, min_length=0)


if __name__ == "__main__":
    unittest.main()
