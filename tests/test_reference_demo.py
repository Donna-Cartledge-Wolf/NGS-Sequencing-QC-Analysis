import tempfile
import unittest
from pathlib import Path
from ngs_qc.reference_demo import reverse_complement, parse_fasta, evaluate_read, run

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data/synthetic_references.fasta"
READS = ROOT / "data/synthetic_alignment_reads.fastq"

class ReferenceDemoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs = parse_fasta(REF)
        from ngs_qc.fastq import read_fastq
        cls.reads = {r.name:r.sequence for r in read_fastq(READS)}

    def call(self, name, **kwargs):
        return evaluate_read(name, self.reads[name], self.refs, **kwargs)

    def test_viral_like_unique(self):
        r=self.call("VIRAL_LIKE_UNIQUE")
        self.assertEqual((r["assignment"],r["best_reference"],r["best_mismatches"]),
                         ("UNIQUE_BEST_IN_TOY_REFERENCES","TOY_VIRAL_LIKE",0))

    def test_host_like_unique(self):
        self.assertEqual(self.call("HOST_LIKE_UNIQUE")["best_reference"],"TOY_HOST_LIKE")

    def test_shared_cross_reference_ambiguity(self):
        r=self.call("AMBIGUOUS_SHARED")
        self.assertEqual(r["assignment"],"CROSS_REFERENCE_AMBIGUOUS")
        self.assertEqual(set(r["equally_best_references"].split(';')),{"TOY_VIRAL_LIKE","TOY_HOST_LIKE"})
        self.assertEqual(r["mapq"],"NOT_CALCULATED")

    def test_reverse_strand(self):
        self.assertEqual(self.call("VIRAL_LIKE_REVERSE")["best_strand"],'-')

    def test_mismatch_identity(self):
        r=self.call("ONE_MISMATCH")
        self.assertEqual((r["best_mismatches"],r["best_identity_pct"]),(1,96.875))

    def test_unambiguous_only(self):
        self.assertEqual(self.call("NON_ACGT")["assignment"],"NON_ACGT")

    def test_short(self):
        self.assertEqual(self.call("TOO_SHORT")["assignment"],"TOO_SHORT")

    def test_no_match(self):
        self.assertEqual(self.call("UNMATCHED")["assignment"],"NO_MATCH_WITHIN_LIMIT")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AGTC"),"GACT")

    def test_bad_threshold(self):
        with self.assertRaises(ValueError):
            evaluate_read('x','ACGT',self.refs,max_mismatches=-1)

    def test_duplicate_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'f.fa'
            p.write_text('>A\nACGT\n>A\nTGCA\n')
            with self.assertRaises(ValueError):parse_fasta(p)

    def test_reports_generated(self):
        with tempfile.TemporaryDirectory() as d:
            s=run(READS,REF,d)
            self.assertEqual(s['reads'],8)
            self.assertTrue((Path(d)/'synthetic_reference_matches.csv').is_file())
