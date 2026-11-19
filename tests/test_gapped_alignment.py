"""Hand-checkable independent expectations for fictional pairwise alignments."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from ngs_qc.gapped_alignment import Scoring, align, read_pairs, run

DATA = Path(__file__).resolve().parents[1] / "data" / "synthetic_gapped_pairs.csv"


class GappedAlignmentTest(unittest.TestCase):
    def test_exact_global(self):
        a = align("ACGT", "ACGT")
        self.assertEqual((a.score, a.matches, a.identity_pct, a.operations), (8, 4, 100.0, "4="))
        self.assertEqual((a.reference_start_1based,a.reference_end_1based), (1,4))

    def test_one_mismatch(self):
        a = align("ACGT", "AGGT")
        self.assertEqual((a.matches, a.substitutions, a.score, a.operations), (3,1,5,"1=1X2="))
        self.assertEqual(a.identity_pct,75.0)

    def test_one_insertion(self):
        a = align("ACGT", "ACGAT")
        self.assertEqual((a.matches,a.insertions_in_query,a.deletions_in_query), (4,1,0))
        self.assertEqual((a.score,a.identity_pct), (6,80.0))
        self.assertIn("-",a.aligned_reference)
        self.assertNotIn("-",a.aligned_query)

    def test_one_deletion(self):
        a = align("ACGT", "AGT")
        self.assertEqual((a.matches,a.insertions_in_query,a.deletions_in_query), (3,0,1))
        self.assertEqual((a.score,a.identity_pct), (4,75.0))
        self.assertNotIn("-",a.aligned_reference)
        self.assertIn("-",a.aligned_query)

    def test_local_shared_motif(self):
        a = align("TTTTACCTGATCGGGG", "GGACCTGATCAAAA", "local")
        self.assertEqual(a.score,16)
        self.assertEqual(a.aligned_reference,"ACCTGATC")
        self.assertEqual(a.aligned_query,"ACCTGATC")
        self.assertEqual((a.reference_start_1based,a.reference_end_1based), (5,12))
        self.assertEqual((a.query_start_1based,a.query_end_1based), (3,10))
        self.assertEqual(a.identity_pct,100.0)

    def test_local_not_global(self):
        ref,query="TTTTACCTGATCGGGG", "GGACCTGATCAAAA"
        local=align(ref,query,"local")
        glob=align(ref,query,"global")
        self.assertGreater(local.score, glob.score)
        self.assertLess(len(local.aligned_query.replace('-','')),len(query))

    def test_local_zero(self):
        a = align("AAAA", "CCCC", "local")
        self.assertEqual((a.score,a.status,a.operations), (0,"NO_POSITIVE_LOCAL_ALIGNMENT",""))
        self.assertIsNone(a.identity_pct)
        self.assertIsNone(a.reference_start_1based)

    def test_global_no_local_hit(self):
        a=align("AAAA","CCCC","global")
        self.assertEqual((a.score,a.substitutions),(-4,4))
        self.assertEqual(a.identity_pct,0.0)

    def test_global_coordinates(self):
        a=align("ACGT","ACGAT","global")
        self.assertEqual((a.reference_start_1based,a.reference_end_1based,a.query_start_1based,a.query_end_1based),(1,4,1,5))

    def test_homopolymer_gap_deterministic_not_proof_of_placement(self):
        a = align("AAAAA","AAAA")
        self.assertEqual((a.matches,a.deletions_in_query,a.score),(4,1,6))
        self.assertEqual(a.identity_pct,80.0)
        self.assertEqual(len(a.aligned_reference),len(a.aligned_query))

    def test_sensitive_to_gap_scoring(self):
        a=align("ACGT", "AGT", scoring=Scoring(match=2,mismatch=-1,gap=-2))
        b=align("ACGT", "AGT", scoring=Scoring(match=2,mismatch=-1,gap=-5))
        self.assertNotEqual(a.score,b.score)

    def test_mismatch_vs_pair_of_gaps_depends_on_scoring(self):
        ref,query="AAA","AAC"
        defaults=align(ref,query,scoring=Scoring(match=2,mismatch=-1,gap=-2))
        changed=align(ref,query,scoring=Scoring(match=2,mismatch=-3,gap=-1))
        self.assertEqual(defaults.operations,"2=1X")
        self.assertEqual((changed.insertions_in_query,changed.deletions_in_query), (1,1))
        self.assertNotEqual(defaults.operations,changed.operations)

    def test_scoring_rejects_bad_signs(self):
        for scores in [(0,-1,-2),(2,1,-2),(2,-1,0),(2,-1,1)]:
            with self.subTest(scores=scores), self.assertRaises(ValueError):
                Scoring(*scores)

    def test_bad_method(self):
        with self.assertRaises(ValueError):align("ACGT","ACGT","semi-global")

    def test_bad_sequences(self):
        for a,b in [("","ACGT"),("ACGN","ACGT"),("ACGT","A-CG"),("A"*501,"ACGT")]:
            with self.subTest(a=a[:8],b=b),self.assertRaises(ValueError):align(a,b)

    def test_one_optimum_not_unique_mapping(self):
        # A repeated run follows the same deterministic tie priority.
        a=align("AAAAA","AAAA")
        self.assertEqual(a,align("AAAAA","AAAA"))
        self.assertEqual(a.deletions_in_query,1)

    def test_pair_csv(self):
        rows=read_pairs(DATA)
        self.assertEqual(len(rows),8)
        self.assertEqual(rows[2]["pair_id"],"INSERTION_IN_QUERY")

    def test_duplicate_pair_rejected(self):
        with TemporaryDirectory() as td:
            p=Path(td)/"dup.csv"
            p.write_text("pair_id,reference,query,teaching_focus\nx,ACGT,ACGT,a\nx,ACGT,ACGT,b\n",encoding="utf8")
            with self.assertRaises(ValueError):read_pairs(p)

    def test_malformed_csv_rejected(self):
        with TemporaryDirectory() as td:
            p=Path(td)/"bad.csv"
            p.write_text("pair_id,reference,query,teaching_focus\nx,ACGT,ACGT\n",encoding="utf8")
            with self.assertRaises(ValueError):read_pairs(p)

    def test_files_and_summary(self):
        with TemporaryDirectory() as td:
            summary=run(DATA, td)
            self.assertEqual((summary["pair_count"],summary["alignment_count"]),(8,16))
            with (Path(td)/"synthetic_gapped_alignments.csv").open(encoding="utf8") as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),16)
            sub=next(row for row in rows if row["pair_id"]=="SUBSTITUTION" and row["method"]=="global")
            self.assertEqual(int(sub["substitutions"]),1)
            d=json.loads((Path(td)/"synthetic_gapped_summary.json").read_text())
            self.assertEqual(d["scoring"]["gap"],-2)


if __name__ == "__main__":
    unittest.main()
