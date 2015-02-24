import tempfile
import unittest
import sqlite3
import os
from StringIO import StringIO

from brocclib.taxonomy_db import (
    init_db, parse_gi_taxid, insert_taxid, parse_names, insert_names,
    parse_nodes, insert_nodes,
    )

class FunctionTests(unittest.TestCase):
    def setUp(self):
        _, self.db = tempfile.mkstemp()
        init_db(self.db)

    def tearDown(self):
        os.remove(self.db)
    
    def test_init_db(self):
        conn = sqlite3.connect(self.db)
        obs = list(conn.execute(
            'SELECT name FROM sqlite_master WHERE type = "table"'))
        self.assertEqual([x[0] for x in obs], ["gi_taxid", "nodes", "names"])

    def test_parse_gi_taxid(self):
        f = StringIO(gi_taxid)
        obs = list(parse_gi_taxid(f))
        exp = [["2", "9913"], ["3", "9913"], ["4", "9646"], ["5", "9913"]]
        self.assertEqual(obs, exp)

    def test_insert_taxid(self):
        f = StringIO(gi_taxid)
        insert_taxid(self.db, f)
        conn = sqlite3.connect(self.db)
        obs = list(conn.execute('SELECT * FROM gi_taxid'))
        exp = [(2, 9913), (3, 9913), (4, 9646), (5, 9913)]
        self.assertEqual(obs, exp)
        
    def test_parse_names(self):
        f = StringIO(names)
        obs = list(parse_names(f))
        self.assertEqual(obs, [["1", "root"], ["2", "Bacteria"]])

    def test_insert_names(self):
        f = StringIO(names)
        insert_names(self.db, f)
        conn = sqlite3.connect(self.db)
        obs = list(conn.execute('SELECT * FROM names'))
        exp = [(1, "root"), (2, "Bacteria")]
        self.assertEqual(obs, exp)

    def test_parse_nodes(self):
        f = StringIO(nodes)
        obs = list(parse_nodes(f))
        self.assertEqual(obs, [
            ["1", "1", "no rank"],
            ["2", "131567", "superkingdom"],
            ["6", "335928", "genus"]])

    def test_insert_nodes(self):
        f = StringIO(nodes)
        insert_nodes(self.db, f)
        conn = sqlite3.connect(self.db)
        obs = list(conn.execute('SELECT * FROM nodes'))
        exp = [
            (1, 1, "no rank"),
            (2, 131567, "superkingdom"),
            (6, 335928, "genus")]
        self.assertEqual(obs, exp)


gi_taxid = """\
2	9913
3	9913
4	9646
5	9913
"""

names = """\
1	|	all	|		|	synonym	|
1	|	root	|		|	scientific name	|
2	|	Bacteria	|	Bacteria <prokaryote>	|	scientific name	|
2	|	Monera	|	Monera <Bacteria>	|	in-part	|
2	|	Procaryotae	|	Procaryotae <Bacteria>	|	in-part	|
2	|	not Bacteria Haeckel 1894	|		|	synonym	|
"""

nodes = """\
1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
2	|	131567	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
6	|	335928	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
"""

if __name__ == "__main__":
    unittest.main()
