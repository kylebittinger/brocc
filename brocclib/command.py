from __future__ import division

import logging
import os

from brocclib.assign import Assigner
from brocclib import get_xml
from brocclib.parse import iter_fasta, read_blast
from brocclib.args import parse_args


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


def main(argv=None):
    opts = parse_args(argv)

    # Configure
    
    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    
    taxa_db = get_xml.NcbiEutils(opts.cache_fp)
    taxa_db.load_cache()

    consensus_thresholds = [t for _, t in CONSENSUS_THRESHOLDS]
    assigner = Assigner(
        opts.min_cover, opts.min_species_id, opts.min_genus_id, opts.min_id,
        consensus_thresholds, opts.max_generic, taxa_db)

    # Read input files
    
    with open(opts.fasta_file) as f:
        sequences = list(iter_fasta(f))

    with open(opts.blast_file) as f:
        blast_hits = read_blast(f)

    # Do the work
        
    assignments = []
    for name, seq in sequences:
        seq_hits = blast_hits[name]
        # This is where the magic happens
        a = assigner.assign(name, seq, seq_hits)
        assignments.append(a)

    # Write output files
        
    taxa_db.save_cache()

    write_output(assignments, opts.output_directory)


def write_output(assignments, output_directory):
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    output_file = open(os.path.join(output_directory, "Full_Taxonomy.txt"), 'w')
    standard_taxa_file = open(os.path.join(output_directory, "Standard_Taxonomy.txt"), "w")
    log_file = open(os.path.join(output_directory, "brocc.log"), "w")
    log_file.write(
        "Sequence\tWinner_Votes\tVotes_Cast\tGenerics_Pruned\tLevel\t"
        "Classification\n")    

    for a in assignments:
        output_file.write(a.format_for_full_taxonomy())
        standard_taxa_file.write(a.format_for_standard_taxonomy())
        log_file.write(a.format_for_log())

    output_file.close()
    standard_taxa_file.close()
    log_file.close()

