import unittest

from brocclib.taxonomy import Lineage

class TaxonTests(unittest.TestCase):
    def setUp(self):
        self.d = [
            ("cellular organisms", "no rank"),
            ("Eukaryota", "superkingdom"),
            ("Opisthokonta", "no rank"),
            ("Fungi", "kingdom"),
            ("Dikarya", "subkingdom"),
            ("Ascomycota", "phylum"),
            ("saccharomyceta", "no rank"),
            ("Saccharomycotina", "subphylum"),
            ("Saccharomycetes", "class"),
            ("Saccharomycetales", "order"),
            # No family
            ("mitosporic Saccharomycetales", "no rank"),
            ("Candida", "genus"),
            ("Candida albicans", "species"),
        ]

    def test_missing_family(self):
        t = Lineage(self.d)
        self.assertEqual(t.get_taxon("family"), "Candida (family)")

    def test_species(self):
        t = Lineage(self.d)
        self.assertEqual(t.get_taxon("species"), "Candida albicans")

    def test_standard_taxa(self):
        t = Lineage(self.d)
        expected = [
            "Eukaryota", "Fungi", "Ascomycota", "Saccharomycetes",
            "Saccharomycetales", "Candida (family)", "Candida",
            "Candida albicans"]
        self.assertEqual(list(t.get_standard_taxa("species")), expected)

    def test_all_taxa(self):
        t = Lineage(self.d)
        expected = [
            "cellular organisms", "Eukaryota", "Opisthokonta",
            "Fungi", "Dikarya", "Ascomycota", "saccharomyceta", "Saccharomycotina",
            "Saccharomycetes", "Saccharomycetales", "mitosporic Saccharomycetales", 
            "Candida", "Candida albicans"]
        self.assertEqual(list(t.get_all_taxa("species")), expected)

    def test_generic_species(self):
        t = Lineage(self.d)
        self.assertEqual(t.is_generic("species"), False)

if __name__ == "__main__":
    unittest.main()

