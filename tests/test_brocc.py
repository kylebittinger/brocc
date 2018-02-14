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
            "--cache_fp", data_fp("brocc_cache.json"),
            ])

    @property
    def _assignments_fp(self):
        return os.path.join(
            self.output_dir, "Standard_Taxonomy.txt")

    def test_no_cache(self):
        main([
            "-i", data_fp("em_10.fasta"),
            "-b", data_fp("em_10_blast.txt"),
            "-o", self.output_dir,
            "-a", "ITS",
        ])
        # A family assignment has changed.
        em_assign = read_from(data_fp("em_10_assignments.txt"))
        em_assign = [a.replace("Coriolaceae", "Steccherinaceae") for a in em_assign]
        self.assertEqual(read_from(self._assignments_fp), em_assign)


    def test_em_10(self):
        self._run_brocc("em_10.fasta", "em_10_blast.txt")
        self.assertEqual(
            read_from(self._assignments_fp),
            read_from(data_fp("em_10_assignments.txt")))

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
