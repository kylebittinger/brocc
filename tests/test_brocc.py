import os.path
import shutil
import tempfile
import unittest

from brocclib.command import main

def data_fp(filename):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'data', filename)


def read_from(filepath):
    with open(filepath) as f:
        res = f.readlines()
    return res


class BroccAcceptance(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp(prefix="brocc")
        
    def tearDown(self):
        shutil.rmtree(self.output_dir)

    def _run_brocc(self, fasta_fp, blast_fp):
        return main([
            "-i", data_fp(fasta_fp),
            "-b", data_fp(blast_fp),
            "-o", self.output_dir,
            "-a", "ITS",
            ])

    @property
    def _assignments_fp(self):
        return os.path.join(
            self.output_dir, "Standard_Taxonomy.txt")

    def test_ize(self):
        self._run_brocc("ize85.fasta", "ize85_blast.txt")
        self.assertEqual(
            read_from(self._assignments_fp),
            read_from(data_fp("ize85_assignments.txt")))

    def test_em_10(self):
        self._run_brocc("em_10.fasta", "em_10_blast.txt")
        self.assertEqual(
            read_from(self._assignments_fp),
            read_from(data_fp("em_10_assignments_2018.txt")))

    def test_em_444(self):
        self._run_brocc("em_444.fasta", "em_444_blast.txt")
        self.assertEqual(
            read_from(self._assignments_fp),
            read_from(data_fp("em_444_assignments.txt")))

    def test_em_312(self):
        self._run_brocc("em_312.fasta", "em_312_blast.txt")
        self.assertEqual(
            read_from(self._assignments_fp),
            read_from(data_fp("em_312_assignments.txt")))

    def test_sac_otu(self):
        self._run_brocc("sac_otu.fasta", "sac_otu_blast.txt")
        self.assertEqual(
            read_from(self._assignments_fp),
            read_from(data_fp("sac_otu_assignments.txt")))

if __name__ == "__main__":
    unittest.main()
