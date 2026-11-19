"""Validation checks for demo arithmetic and defensive input handling."""
import csv
import tempfile
import unittest
from pathlib import Path
from ngs_qc.control_review import load_records, make_review, run

FIXTURE = Path(__file__).resolve().parents[1] / 'data' / 'synthetic_control_summary.csv'

class ReviewTest(unittest.TestCase):
    def setUp(self):
        self.records = load_records(FIXTURE)
        self.rows = make_review(self.records)

    def test_four_sample_target_rows(self):
        self.assertEqual(len(self.rows), 4)

    def test_sample_a_rpm_and_ratio(self):
        item = next(r for r in self.rows if r['sample_id'] == 'CLIENT_A')
        self.assertAlmostEqual(item['sample_rpm'], 5.833333, places=5)
        self.assertEqual(item['blank_rpm'], 4.0)
        self.assertAlmostEqual(item['sample_to_blank_rpm_ratio'], 1.458333, places=5)

    def test_blank_zero_leads_to_undefined_ratio(self):
        item = next(r for r in self.rows if r['sample_id'] == 'CLIENT_C')
        self.assertIsNone(item['sample_to_blank_rpm_ratio'])

    def test_zero_reads_never_called_negative(self):
        item = next(r for r in self.rows if r['sample_id'] == 'CLIENT_B')
        self.assertEqual(item['sample_putative_reads'], 0)
        self.assertEqual(item['review_status'], 'TECHNICAL_REVIEW_REQUIRED')
        self.assertIn('does not demonstrate virus absence', item['review_notes'])

    def test_no_positive_or_negative_classification(self):
        self.assertTrue(all(r['review_status'] == 'TECHNICAL_REVIEW_REQUIRED' for r in self.rows))

    def test_csv_and_json_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = run(FIXTURE, tmp)
            self.assertEqual(len(out['results']), 4)
            self.assertTrue((Path(tmp) / 'viral_control_review.csv').exists())
            self.assertTrue((Path(tmp) / 'viral_control_review.json').exists())

    def test_missing_blank_rejected(self):
        no_blank = [r for r in self.records if r.sample_type != 'extraction_blank']
        with self.assertRaisesRegex(ValueError, 'No extraction blank'):
            make_review(no_blank)

    def test_duplicate_blank_rejected(self):
        blank = next(r for r in self.records if r.sample_type == 'extraction_blank')
        with self.assertRaisesRegex(ValueError, 'Duplicate extraction blanks'):
            make_review(self.records + [blank])

    def test_duplicate_sample_row_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / 'duplicate.csv'
            text = FIXTURE.read_text()
            dest.write_text(text + text.splitlines()[1] + '\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                load_records(dest)

    def test_zero_depth_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / 'zero.csv'
            lines = FIXTURE.read_text().splitlines()
            values = next(csv.reader([lines[1]]))
            values[4] = '0'
            lines[1] = ','.join(values)
            dest.write_text('\n'.join(lines) + '\n')
            with self.assertRaisesRegex(ValueError, 'impossible'):
                load_records(dest)

if __name__ == '__main__':
    unittest.main()
