import unittest

from brocclib.taxonomy import Lineage

class TaxonTests(unittest.TestCase):
    def setUp(self):
        self.d = {
            "species": "Candida albicans",
            "genus": "Candida",
            # no family
            "order": "Saccharomycetales",
            "class": "Saccharomycetes",
            "phylum": "Ascomycota",
            "kingdom": "Fungi",
            "superkingdom": "Eukaryota",
            "Lineage": (
                "cellular organisms; Eukaryota; Opisthokonta; Fungi; "
                "Dikarya; Ascomycota; saccharomyceta; Saccharomycotina; "
                "Saccharomycetes; Saccharomycetales; mitosporic "
                "Saccharomycetales; Candida"),
            "LineageWithRanks": [
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
                ("mitosporic Saccharomycetales", "no rank"),
                ("Candida", "genus"),
                ("Candida albicans", "species"),
            ]
        }

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
        self.assertEqual(t.classified, True)

    def test_all_taxa(self):
        t = Lineage(self.d)
        expected = [
            "cellular organisms", "Eukaryota", "Opisthokonta",
            "Fungi", "Dikarya", "Ascomycota", "saccharomyceta", "Saccharomycotina",
            "Saccharomycetes", "Saccharomycetales", "mitosporic Saccharomycetales", 
            "Candida", "Candida albicans"]
        self.assertEqual(list(t.get_all_taxa("species")), expected)
        self.assertEqual(t.classified, True)

    def test_generic_species(self):
        self.d["species"] = "uncultured organism"
        t = Lineage(self.d)
        self.assertEqual(t.classified, False)

    def test_generic_taxon(self):
        self.d["no rank"] = "unclassified Fungi"
        t = Lineage(self.d)
        self.assertEqual(t.classified, False)


if __name__ == "__main__":
    unittest.main()

