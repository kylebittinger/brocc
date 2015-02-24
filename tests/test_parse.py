from unittest import TestCase, main
from cStringIO import StringIO

from brocclib.parse import (
    read_blast, iter_fasta, parse_gi_number,
    )


class GiNumberTests(TestCase):
    def test_gi_start(self):
        obs = parse_gi_number('gi|259100874|gb|GQ513762.1|')
        self.assertEqual(obs, '259100874')

    def test_gi_empty(self):
        for s in ['', 'gi', 'ran|dom']:
            self.assertEqual(parse_gi_number(s), None)


class FastaTests(TestCase):
    def test_basic(self):
        lines = [
                 ">lab1",
                 "TTTTCCC",
                 ">lab2",
                 "CCAAAA",
                 ]
        seqs = iter_fasta(lines)
        self.assertEqual(seqs.next(), ("lab1", "TTTTCCC"))
        self.assertEqual(seqs.next(), ("lab2", "CCAAAA"))
        self.assertRaises(StopIteration, seqs.next)


class BlastOutputTests(TestCase):
    def test_normal_output(self):
        obs = read_blast(StringIO(normal_output))
        h = obs['0 E7_168192'][0]
        self.assertEqual(h.gi, "259100874")
        self.assertEqual(h.pct_id, 98.74)
        self.assertEqual(h.length, 159)

    def test_malformed_output(self):
        obs = read_blast(StringIO(malformed_output))
        h = obs['0 E7_168192'][0]
        self.assertEqual(h.gi, "259100874")
        self.assertEqual(h.pct_id, 98.74)
        self.assertEqual(h.length, 159)

    def test_missing_read(self):
        obs = read_blast(StringIO(normal_output))
        self.assertEqual(obs['sdlkj'], [])

    


normal_output = """\
# BLASTN 2.2.25+
# Query:  0 E7_168192
# Database: /home/rohinis/blastdb/blast_nt/nt
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 100 hits found
0	gi|259100874|gb|GQ513762.1|	98.74	159	1	1	407	564	1	159	2e-70	 275
0	gi|259098555|gb|GQ520853.1|	98.74	159	1	1	407	564	1	159	2e-70	 275
0	gi|259098210|gb|GQ520508.1|	98.11	159	2	1	407	564	1	159	1e-68	 269
0	gi|259092808|gb|GQ524514.1|	98.11	159	1	2	407	564	1	158	1e-67	 266
0	gi|259107208|gb|GQ510686.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
0	gi|259103360|gb|GQ516248.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
0	gi|259101730|gb|GQ514618.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
# BLAST processed 608 queries
"""


malformed_output = """\
# BLASTN 2.2.25+
# Query: 0 E7_168192
# Database: /home/rohinis/blastdb/blast_nt/nt
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 100 hits found
0	gi|259100874|gb|GQ513762.1|	98.74	159	1	1	407	564	1	159	2e-70	 275
0	gi|259098555|gb|GQ520853.1|	98.74	159	1	1	407	564	1	159	2e-70	 275
0	gi|259098210|gb|GQ520508.1|	98.11	159	2	1	407	564	1	159	1e-68	 269
0	gi|259092808|gb|GQ524514.1|	98.11	159	1	2	407	564	1	158	1e-67	 266
0	gi|259107208|gb|GQ510686.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
0	gi|259103360|gb|GQ516248.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
0	gi|259101730|gb|GQ514618.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
0	gi|259093119|gb|GQ524825.1|	98.68	152	1	1	414	564	1	152	2e-66	 262
0	gi|259100068|gb|GQ522366.1|	98.67	150	1	1	416	564	1	150	2e-65	 259
0	gi|259099396|gb|GQ521694.1|	98.67	150	1	1	416	564	1	150	2e-65	 259
"""


if __name__ == '__main__':
    main()
