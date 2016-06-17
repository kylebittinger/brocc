import tempfile
import unittest

from brocclib.get_xml import (
    get_taxid, get_lineage, NcbiEutils, get_taxid_from_accession,
    )


class NcbiEutilsTests(unittest.TestCase):
    def setUp(self):
        self.cache_file = tempfile.NamedTemporaryFile(suffix=".json")
        self.db = NcbiEutils(self.cache_file.name)

    def test_save_load_cache(self):
        lineages = {
            "taxon1": {'class': "a", "genus": "b"},
            "taxon2": {'class': "c", "genus": "d"},
            }
        taxon_ids = {"taxon1": "b", "taxon2": "d"}
        self.db.lineages = lineages
        self.db.taxon_ids = taxon_ids
        self.db._fresh = False
        self.db.save_cache()

        db2 = NcbiEutils(self.cache_file.name)
        db2.load_cache()
        self.assertEqual(db2.lineages, lineages)
        self.assertEqual(db2.taxon_ids, taxon_ids)

    def test_get_taxon_id(self):
        self.assertEqual(self.db.get_taxon_id("312434489"), "531911")
        self.assertEqual(self.db.taxon_ids, {"312434489": "531911"})

    def test_get_lineage(self):
        observed_lineage = self.db.get_lineage("531911")
        expected_lineage = {
            'superkingdom': 'Eukaryota',
            'kingdom': 'Fungi',
            'no rank': 'sordariomyceta',
            'family': 'Sporocadaceae',
            'Lineage': (
                'cellular organisms; Eukaryota; Opisthokonta; Fungi; Dikarya; '
                'Ascomycota; saccharomyceta; Pezizomycotina; leotiomyceta; '
                'sordariomyceta; Sordariomycetes; Xylariomycetidae; '
                'Xylariales; Sporocadaceae; Pestalotiopsis'),
            'subkingdom': 'Dikarya',
            'subclass': 'Xylariomycetidae',
            'order': 'Xylariales',
            'phylum': 'Ascomycota',
            'species': 'Pestalotiopsis maculiformans',
            'subphylum': 'Pezizomycotina',
            'genus': 'Pestalotiopsis',
            'class': 'Sordariomycetes'}
        self.assertEqual(observed_lineage, expected_lineage)
        self.assertEqual(self.db.lineages, {'531911': expected_lineage})


class FunctionTests(unittest.TestCase):
    def test_get_taxid(self):
        self.assertEqual(get_taxid("312434489"), "531911")

    def test_get_taxid_from_accession(self):
        self.assertEqual(get_taxid("HQ844023.1"), "1056490")

    def test_getLineage(self):
        lineage = get_lineage("531911")
        self.assertEqual(lineage["genus"], "Pestalotiopsis")
        self.assertEqual(lineage["kingdom"], "Fungi")
        # Should this return the HTTP 400 error?
        self.assertEqual(get_lineage("asdf"), None)


if __name__ == '__main__':
    unittest.main()
