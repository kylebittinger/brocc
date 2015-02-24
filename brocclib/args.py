import optparse

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
    parser.add_option("--max_generic", type="float", default=.7, help=(
        "maximum proportion of generic classifications allowed "
        "before query cannot be classified [default: %default]"))
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
            
