import optparse
import sqlite3
import os.path
import gzip
import itertools
import tarfile

# Columns kept in node table: 1, 2, 3

SCHEMA = """\
CREATE TABLE gi_taxid(
       nuc_id INT NOT NULL PRIMARY KEY,
       tax_id INT NOT NULL);

CREATE TABLE nodes (
  tax_id INT NOT NULL PRIMARY KEY,
  parent_id INT NOT NULL,
  rank TEXT);

CREATE TABLE names (
  tax_id INT NOT NULL PRIMARY KEY,
  name TEXT);
"""

def init_db(fp):
    """Create a new SQLite3 database for the NCBI taxonomy.
    """
    conn = sqlite3.connect(fp)
    conn.executescript(SCHEMA)
    conn.commit()
    conn.close()


def parse_gi_taxid(f):
    for line in f:
        line = line.rstrip()
        if line:
            yield line.split("\t")


def _parse_ncbi_dmp(f):
    for line in f:
        line = line.rstrip("\t|\n")
        if line:
            yield line.split("\t|\t")


def parse_names(f):
    rows = _parse_ncbi_dmp(f)
    return (r[0:2] for r in rows if r[3] == "scientific name")


def parse_nodes(f):
    rows = _parse_ncbi_dmp(f)
    return (r[0:3] for r in rows)


def _insert_many(db, input_file, parse_fcn, insert_sql):
    rows = parse_fcn(input_file)
    conn = sqlite3.connect(db)
    conn.executemany(insert_sql, rows)
    conn.commit()
    conn.close()


def insert_names(db, f):
    sql = "INSERT INTO names VALUES (?,?)"
    return _insert_many(db, f, parse_names, sql)


def insert_nodes(db, f):
    sql = "INSERT INTO nodes VALUES (?,?,?)"
    return _insert_many(db, f, parse_nodes, sql)


def insert_taxid(db, f):
    sql = "INSERT INTO gi_taxid VALUES (?,?)"
    return _insert_many(db, f, parse_gi_taxid, sql)


def main(argv=None):
    p = optparse.OptionParser()
    p.add_option("--taxid_fp", help="Path to gzipped taxid file")
    p.add_option("--taxdmp_fp", help="Path to tar-gzipped taxdmp file")
    p.add_option("--db_fp", help="Output filepath for sqlite3 database")
    opts, args = p.parse_args(argv)

    if os.path.exists(opts.db_fp):
        p.error("Database file already exists.  Please delete first.")
    if not os.path.exists(opts.taxid_fp):
        p.error("Taxid file not found.")
    if not os.path.exists(opts.taxdmp_fp):
        p.error("Taxdmp file not found.")

    taxid_f = gzip.GzipFile(opts.taxid_fp)
    taxdmp_tar = tarfile.open(opts.taxdmp_fp)
    nodes_f = taxdmp_tar.extractfile("nodes.dmp")
    names_f = taxdmp_tar.extractfile("names.dmp")

    init_db(opts.db_fp)
    insert_taxid(opts.db_fp, taxid_f)
    insert_nodes(opts.db_fp, nodes_f)
    insert_names(opts.db_fp, names_f)
