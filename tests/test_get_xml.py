import tempfile
import unittest

from brocclib.get_xml import (
    get_taxid, get_lineage, NcbiEutils,
    )


class NcbiEutilsTests(unittest.TestCase):
    def test_get_taxon_id(self):
        db = NcbiEutils()
        self.assertEqual(db.get_taxon_id("HQ608011.1"), "531911")
        self.assertEqual(db.taxon_ids, {"HQ608011.1": "531911"})

    def test_get_lineage(self):
        db = NcbiEutils()
        observed_lineage = db.get_lineage("531911")
        expected_lineage = [
                ('cellular organisms', 'no rank'),
                ('Eukaryota', 'superkingdom'),
                ('Opisthokonta', 'no rank'),
                ('Fungi', 'kingdom'),
                ('Dikarya', 'subkingdom'),
                ('Ascomycota', 'phylum'),
                ('saccharomyceta', 'no rank'),
                ('Pezizomycotina', 'subphylum'),
                ('leotiomyceta', 'no rank'),
                ('sordariomyceta', 'no rank'),
                ('Sordariomycetes', 'class'),
                ('Xylariomycetidae', 'subclass'),
                ('Xylariales', 'order'),
                ('Sporocadaceae', 'family'),
                ('Pestalotiopsis', 'genus'),
                ('Pestalotiopsis maculiformans', 'species')]
        self.assertEqual(observed_lineage, expected_lineage)
        self.assertEqual(db.lineages, {'531911': expected_lineage})


class FunctionTests(unittest.TestCase):
    def test_get_taxid(self):
        self.assertEqual(get_taxid("312434489"), "531911")

    def test_get_taxid_from_accession(self):
        self.assertEqual(get_taxid("HQ844023.1"), "1056490")

    def test_getLineage(self):
        # Should this return the HTTP 400 error?
        self.assertEqual(get_lineage("asdf"), None)


if __name__ == '__main__':
    unittest.main()
