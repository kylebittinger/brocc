from __future__ import division

import logging
import optparse
import os
import shutil
import subprocess
import tempfile

from brocclib.assign import Assigner
from brocclib.get_xml import NcbiEutils
from brocclib.parse import iter_fasta, read_blast


'''
Created on Aug 29, 2011
@author: Serena, Kyle
'''


CONSENSUS_THRESHOLDS = [
    ("species", 0.6),
    ("genus", 0.6),
    ("family", 0.6),
    ("order", 0.9),
    ("clas", 0.9),
    ("phylum", 0.9),
    ("kingdom", 0.9),
    ("domain", 0.9),
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
    parser.add_option("--min_votes", type="int", default=10, help=(
        "minimum number of votes needed to establish a consensus "
        "after removal of generic taxa [default: %default]"))
    parser.add_option("--cache_fp", help=(
        "Filepath for retaining data retrieved from NCBI between runs.  "
        "Can help to reduce execution time if BROCC is run several times."))
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
    
    taxa_db = NcbiEutils(opts.cache_fp)
    taxa_db.load_cache()

    consensus_thresholds = [t for _, t in CONSENSUS_THRESHOLDS]
    assigner = Assigner(
        opts.min_cover, opts.min_species_id, opts.min_genus_id, opts.min_id,
        consensus_thresholds, opts.min_votes, taxa_db)

    # Read input files
    
    with open(opts.fasta_file) as f:
        sequences = list(iter_fasta(f))

    with open(opts.blast_file) as f:
        blast_hits = read_blast(f)

    # Open output files

    if not os.path.exists(opts.output_directory):
        os.mkdir(opts.output_directory)
    output_file = open(
        os.path.join(opts.output_directory, "Full_Taxonomy.txt"), 'w')
    standard_taxa_file = open(
        os.path.join(opts.output_directory, "Standard_Taxonomy.txt"), "w")
    log_file = open(os.path.join(opts.output_directory, "brocc.log"), "w")
    log_file.write(
        "Sequence\tWinner_Votes\tVotes_Cast\tGenerics_Pruned\tLevel\t"
        "Classification\n")

    # Do the work

    for name, seq in sequences:
        seq_hits = blast_hits[name]
        # This is where the magic happens
        a = assigner.assign(name, seq, seq_hits)

        output_file.write(a.format_for_full_taxonomy())
        standard_taxa_file.write(a.format_for_standard_taxonomy())
        log_file.write(a.format_for_log())

    # Close output files, write cache

    output_file.close()
    standard_taxa_file.close()
    log_file.close()

    taxa_db.save_cache()

def run_acceptance_test(argv=None):
    p = optparse.OptionParser()
    p.add_option("--no_cache", action="store_true")
    opts, args = p.parse_args(argv)

    base_file_paths = [os.path.splitext(fp)[0] for fp in args]
    for base_fp in set(base_file_paths):
        fasta_fp = "{0}.fasta".format(base_fp)
        blast_fp = "{0}_blast.txt".format(base_fp)

        output_dir = tempfile.mkdtemp(prefix="brocc")

        print "Temporary output directory:", output_dir

        brocc_args = [
            "-i", fasta_fp, "-b", blast_fp, "-o", output_dir,
            "-a" "ITS", "--verbose"]
        if not opts.no_cache:
            cache_fp = os.path.expanduser("~/.brocc_acceptance_tests.json")
            print "Using cache", cache_fp
            brocc_args.extend(["--cache_fp", cache_fp])
        main(brocc_args)

        observed_assignments_fp = os.path.join(
            output_dir, "Standard_Taxonomy.txt")
        expected_assignments_fp = "{0}_assignments.txt".format(base_fp)
        subprocess.call(
            ["diff", observed_assignments_fp, expected_assignments_fp],
        )
        #shutil.rmtree(output_dir)
