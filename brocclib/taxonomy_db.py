import argparse
import gzip
import itertools
import optparse
import os
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile

ACCESSION_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
TAXDUMP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

CONFIG_DIR = os.path.expanduser("~/.brocc")
TAXONOMY_DB_FILENAME = "taxonomy.db"
TAXONOMY_DB_FP = os.path.join(CONFIG_DIR, TAXONOMY_DB_FILENAME)

def parse_names(f):
    for rec in parse_ncbi_table(f):
        taxid, name, _, name_class = rec
        if name_class == "scientific name":
            yield taxid, name

def parse_nodes(f):
    for rec in parse_ncbi_table(f):
        taxid = rec[0]
        parent = rec[1]
        rank = rec[2]
        yield taxid, parent, rank

def parse_ncbi_table(f):
    for line in f:
        yield line.rstrip("\t|\n").split("\t|\t")

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--download_dir",
        help="Directory to download files from NCBI. (default: temp dir)")
    p.add_argument(
        "--sqlite_db", default=TAXONOMY_DB_FP,
        help="Path to sqlite3 database (default: %(default)s)")
    args = p.parse_args(argv)

    if args.download_dir is None:
        download_dir = tempfile.mkdtemp()
    else:
        download_dir = args.download_dir
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)

    accession_fp = download_accessions(download_dir)
    names_fp, nodes_fp = download_nodes(download_dir)

    with open(accession_fp) as f_in:
        accessions = list(parse_accessions(f_in))
    with open(names_fp) as f_names, open(nodes_fp) as f_nodes:
        names = parse_names(f_names)
        unnamed_nodes = parse_nodes(f_nodes)
        nodes = list(combine_names_and_nodes(names, unnamed_nodes))

    sys.stderr.write("Initializing sqlite3 database...")
    sqlite_db = os.path.expanduser(args.sqlite_db)
    if os.path.exists(sqlite_db):
        os.remove(sqlite_db)
    database_dir = os.path.dirname(sqlite_db)
    database_dir_is_default = database_dir == CONFIG_DIR
    if database_dir_is_default and (not os.path.exists(CONFIG_DIR)):
        os.mkdir(CONFIG_DIR)
    init_db(sqlite_db, download_dir, accessions, nodes)

    shutil.rmtree(download_dir)

def download_accessions(download_dir):
    zipped_accession_fp = os.path.join(
        download_dir, os.path.basename(ACCESSION_URL))
    accession_fp = os.path.splitext(zipped_accession_fp)[0]
    if not os.path.exists(accession_fp):
        subprocess.check_call(
            ["wget", "--directory-prefix", download_dir, ACCESSION_URL])
        subprocess.check_call(["gunzip", zipped_accession_fp])
    return accession_fp

def download_nodes(download_dir):
    taxdump_fp = os.path.join(download_dir, os.path.basename(TAXDUMP_URL))
    nodes_fp = os.path.join(download_dir, "nodes.dmp")
    names_fp = os.path.join(download_dir, "names.dmp")
    if not all(os.path.exists(x) for x in (nodes_fp, names_fp)):
        if not os.path.exists(taxdump_fp):
            subprocess.check_call(
                ["wget", "--directory-prefix", download_dir, TAXDUMP_URL])
        subprocess.check_call(["tar", "xvzf", taxdump_fp, "-C", download_dir])
    return names_fp, nodes_fp


def combine_names_and_nodes(names, nodes):
    names = dict(names)
    for taxid, parent, rank in nodes:
        name = names[taxid]
        yield taxid, parent, name, rank

def parse_accessions(f):
    # First line is the header
    next(f)
    for line in f:
        vals = line.rstrip().split("\t")
        unversioned_accession = vals[0]
        taxid = vals[2]
        yield unversioned_accession, taxid

def unversion(acc):
    if "." in acc:
        return acc.rpartition(".")[0]
    else:
        return acc

class NcbiLocal(object):
    select_taxon_id = "SELECT taxid FROM accessions WHERE accession = ?"
    select_node = "SELECT parent, name, rank FROM nodes WHERE taxid = ?"

    def __init__(self, db):
        self.db = db
        self.con = sqlite3.connect(self.db)

    def get_taxon_id(self, acc):
        unversioned_acc = unversion(acc)
        cur = self.con.cursor()
        cur.execute(self.select_taxon_id, (unversioned_acc,))
        self.con.commit()
        res = cur.fetchone()
        cur.close()
        if res is None:
            return None
        else:
            return res[0]

    def get_lineage(self, taxon_id):
        max_iter = 100
        lineage = list(self._query_nodes(taxon_id))
        res = {}
        for name, rank in reversed(lineage):
            res[rank] = name
        res["Lineage"] = "; ".join(name for name, rank in reversed(lineage))
        return res

    def _query_nodes(self, taxon_id):
        max_iter = 100
        for num_iter in range(max_iter):
            cur = self.con.cursor()
            cur.execute(self.select_node, (taxon_id,))
            self.con.commit()
            res = cur.fetchone()
            cur.close()
            if res is None:
                break
            parent, name, rank = res
            yield name, rank
            if parent == taxon_id:
                break
            else:
                taxon_id = parent

    def save_cache(self):
        pass

    def load_cache(self):
        pass


def init_db(db, temp_dir, accessions, nodes):
    accession_fp = os.path.join(temp_dir, "accessions_import.txt")
    with open(accession_fp, "w") as f:
        _write_tsv(accessions, f)

    nodes_fp = os.path.join(temp_dir, "nodes_import.txt")
    with open(nodes_fp, "w") as f:
        _write_tsv(nodes, f)

    commands_fp = os.path.join(temp_dir, "create_taxondb.sql")
    with open(commands_fp, "w") as f:
        f.write(SQLITE3_IMPORT_COMMANDS)

    with open(commands_fp) as f:
        subprocess.check_call(["sqlite3", db], stdin=f, cwd=temp_dir)

SQLITE3_IMPORT_COMMANDS = """\
CREATE TABLE accessions (
    "accession" TEXT,
    "taxid" INTEGER
);
CREATE TABLE nodes (
    "taxid" INTEGER,
    "parent" INTEGER,
    "name" TEXT,
    "rank" TEXT
);
.separator "\t"
.import accessions_import.txt accessions
.import nodes_import.txt nodes
CREATE UNIQUE INDEX idx_accessions ON accessions(accession);
CREATE UNIQUE INDEX idx_taxid ON nodes(taxid);
"""


def _write_tsv(recs, f):
    for rec in recs:
        f.write("\t".join(str(r) for r in rec))
        f.write("\n")
