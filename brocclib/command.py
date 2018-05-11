from __future__ import division

import logging
import optparse
import os
import shutil
import subprocess
import sys
import tempfile

from brocclib.assign import Assigner
from brocclib.get_xml import NcbiEutils
from brocclib.taxonomy_db import NcbiLocal, TAXONOMY_DB_FP
from brocclib.parse import iter_fasta, read_blast


'''
Created on Aug 29, 2011
@author: Serena, Kyle
'''


CONSENSUS_THRESHOLDS = [
    ("species", 0.6),
    ("genus", 0.6),
    ("family", 0.6),
    ("order", 0.7),
    ("class", 0.8),
    ("phylum", 0.9),
    ("kingdom", 0.9),
    ("superkingdom", 0.9),
    ]


def parse_args(argv=None):
    parser = optparse.OptionParser(description=(
        "BROCC uses a consensus method determine taxonomic assignments from "
        "BLAST hits."))
    parser.add_option("--min_id", type="float", default=80.0, help=(
        "minimum identity required for a db hit to be considered at any "
        "level [default: %default]"))
    parser.add_option("--min_cover", type="float", default=.7, help=(
        "minimum coverage required for a db hit to be considered "
        "[default: %default]"))
    parser.add_option("--min_species_id", type="float", help=(
        "minimum identity required for a db hit to be "
        "considered at species level [default: %default]"))
    parser.add_option("--min_genus_id", type="float", help=(
        "minimum identity required for a db hit to be "
        "considered at genus level [default: %default]"))
    parser.add_option("--min_winning_votes", type="int", default=4, help=(
        "minimum number of votes needed to establish a consensus "
        "after removal of generic taxa [default: %default]"))
    parser.add_option("--taxonomy_db", default=TAXONOMY_DB_FP, help=(
        "location of sqlite3 database holding a local copy of the "
        "NCBI taxonomy [default: %default]"))
    parser.add_option("-v", "--verbose", action="store_true",
        help="output message after every query sequence is classified")
    parser.add_option("-i", "--input_fasta_file", dest="fasta_file",
        help="input fasta file of query sequences [REQUIRED]")
    parser.add_option("-b", "--input_blast_file", dest="blast_file",
        help="input blast file [REQUIRED]")
    parser.add_option("-o", "--output_directory",
        help="output directory [REQUIRED]")
    parser.add_option("-a", "--amplicon", help=(
        "amplicon being classified, either 'ITS' or '18S'. If this option is "
        "not supplied, both --min_species_id and --min_genus_id must be "
        "specified"))
    opts, args = parser.parse_args(argv)

    if opts.amplicon == "ITS":
        opts.min_genus_id = 83.05
        opts.min_species_id = 95.2
    elif opts.amplicon == "18S":
        opts.min_genus_id = 96.0
        opts.min_species_id = 99.0
    elif opts.amplicon:
        parser.error("Provided amplicon %s not recognized." % opts.amplicon)
    else:
        if not (opts.min_species_id and opts.min_genus_id):
            parser.error("Must specify --amplicon, or provide both --min_species_id and --min_genus_id.")

    return opts


def main(argv=None):
    opts = parse_args(argv)

    # Configure
    
    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    if os.path.exists(opts.taxonomy_db):
        taxa_db = NcbiLocal(opts.taxonomy_db)
    else:
        sys.stderr.write(
            "Did not detect a local copy of the NCBI taxonomy.\n"
            "Using NCBI EUtils to get taxonomic info instead.\n\n"
            "The NCBI taxonomy can be dowloaded with the script "
            "create_local_taxonomy_db.py\n"
            "This will greatly speed up the assignment process.\n"
        )
        taxa_db = NcbiEutils()

    consensus_thresholds = [t for _, t in CONSENSUS_THRESHOLDS]
    assigner = Assigner(
        opts.min_cover, opts.min_species_id, opts.min_genus_id, opts.min_id,
        consensus_thresholds, opts.min_winning_votes, taxa_db)

    # Read input files

    with open(opts.fasta_file) as f:
        sequences = list(iter_fasta(f))

    with open(opts.blast_file) as f:
        blast_hits = read_blast(f)

    # Open output files

    if not os.path.exists(opts.output_directory):
        os.mkdir(opts.output_directory)
    standard_taxa_file = open(
        os.path.join(opts.output_directory, "Standard_Taxonomy.txt"), "w")
    log_file = open(os.path.join(opts.output_directory, "brocc.log"), "w")
    log_file.write(
        "Sequence\tWinner_Votes\tVotes_Cast\tGenerics_Pruned\tLevel\t"
        "Classification\n")

    # Set up log for voting details
    vote_logger = logging.getLogger("brocc.votes")
    vote_logger.setLevel(logging.DEBUG)
    vote_handler = logging.FileHandler(os.path.join(opts.output_directory, "voting_log.txt"))
    vote_handler.setLevel(logging.DEBUG)
    vote_formatter = logging.Formatter('%(message)s')
    vote_handler.setFormatter(vote_formatter)
    vote_logger.addHandler(vote_handler)
    vote_logger.propagate = False
    # Do the work

    for name, seq in sequences:
        seq_hits = blast_hits[name]
        # This is where the magic happens
        a = assigner.assign(name, seq, seq_hits)

        standard_taxa_file.write(a.format_for_standard_taxonomy())
        log_file.write(a.format_for_log())

    # Close output files

    standard_taxa_file.close()
    log_file.close()

def run_comparison(argv=None):
    p = optparse.OptionParser()
    p.add_option("--keep_temp", action="store_true")
    opts, args = p.parse_args(argv)

    base_file_paths = [os.path.splitext(fp)[0] for fp in args]
    for base_fp in set(base_file_paths):
        fasta_fp = "{0}.fasta".format(base_fp)
        blast_fp = "{0}_blast.txt".format(base_fp)

        output_dir = tempfile.mkdtemp(prefix="brocc")

        brocc_args = [
            "-i", fasta_fp, "-b", blast_fp, "-o", output_dir,
            "-a" "ITS"]
        main(brocc_args)

        base_filename = os.path.basename(base_fp)

        voting_src = os.path.join(output_dir, "voting_log.txt")
        voting_dest = "{0}_voting_log.txt".format(base_filename)
        shutil.copyfile(voting_src, voting_dest)

        observed_assignments_fp = os.path.join(
            output_dir, "Standard_Taxonomy.txt")
        expected_assignments_fp = "{0}_assignments.txt".format(base_fp)
        diff_fp = "{0}_diff.txt".format(base_filename)
        with open(diff_fp, "w") as f:
            subprocess.call(
                ["diff", observed_assignments_fp, expected_assignments_fp],
                stdout=f,
            )

        if not opts.keep_temp:
            shutil.rmtree(output_dir)
