import os
import shutil
import tempfile
import unittest

from brocclib.taxonomy_db import (
    NcbiLocal, init_db,
)

TEST_ACCESSIONS = [
    ("ABC123", 3324),
    ("AF56.1", 9012),
]

TEST_NODES = [
    (3324, 5, "Clostridiales", "order"),
    (5, 4, "Clostridia", "class"),
    (4, 3, "Firmicutes", "phylum"),
    (3, 2, "Bacteria", "superkingdom"),
    (2, 335, "cellular organisms", "no rank"),
    (1, 1, "root", "no rank"),
    (9012, 5, "Natronoanaerobium aggerbacterium", "species"),
]

class NcbiLocalTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        sqlite_fp = os.path.join(self.temp_dir, "taxonomy.db")
        init_db(sqlite_fp, self.temp_dir, TEST_ACCESSIONS, TEST_NODES)
        self.db = NcbiLocal(sqlite_fp)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_lineage(self):
        observed_lineage = self.db.get_lineage(3324)
        expected_lineage = {
            u'class': u'Clostridia',
            u'no rank': u'cellular organisms',
            u'order': u'Clostridiales',
            u'phylum': u'Firmicutes',
            u'superkingdom': u'Bacteria',
            'Lineage': (
                u'cellular organisms; Bacteria; Firmicutes; '
                u'Clostridia; Clostridiales'),
        }
        self.assertEqual(observed_lineage, expected_lineage)
