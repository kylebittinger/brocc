BROCC
=====

Generate consensus-based taxonomic assignments from BLAST results.

Please cite the following paper when using BROCC:

Dollive S, Peterfreund GL, Sherrill-Mix S, Bittinger K, Sinha R, Hoffmann C, Nabel CS, Hill DA, Artis D, Bachman MA, Custers-Allen R, Grunberg S, Wu GD, Lewis JD, Bushman FD.  A tool kit for quantifying eukaryotic rRNA gene sequences from human microbiome samples.  Genome Biol. 2012 Jul 3;13(7):R60. doi: 10.1186/gb-2012-13-7-r60.

Installion
----------

To install BROCC, run this command in the current directory.

    pip install .

BROCC has two methods of looking up taxon names: it can use a local
copy of the NCBI taxonomy, or it can use NCBI's E-utilities to ask for
names over the web.  Using a local copy of the NCBI taxonomy is much
faster and more reliable.  To prepare this taxonomy database, use the
command:

    create_local_taxonomy_db.py

You will need about 5G for the taxonomy databaase, which is stored at
`~/.brocc/taxonomy.db` by default.

Running
-------

The BROCC classifier takes BLAST results as input, using output format
7 (see BLAST documentation).  The following BLAST parameters are
found to work best for amplicon-based sequence sets:

    blastn -query <SEQUENCES (FASTA FORMAT)> -evalue 1e-5 -outfmt 7 -db nt -out <BLAST RESULTS> -num_threads 8 -max_target_seqs 100

The BROCC program requires two input files and the name of an output directory:

    brocc.py -i <SEQUENCES (FASTA FORMAT)> -b <BLAST RESULTS> -o <OUTPUT DIRECTORY>

`brocc.py` outputs a QIIME-formated taxonomy map and a log file.  The
log file that contains the full classification and voting details:
number of votes for winner, total votes cast, and number of generic
hits pruned.

Settings
--------

The BROCC command has several options for consensus formation:

* minimum hit coverage (for consideration)
* minimum species identity (for consideration at the species level)
* minimum genus identity (for consideration at the genus level)

The defaults are currently set for the ITS1 gene, because these
settings seem to work well over several different amplicons.  The
minimum identity defaults for ITS1 are 95.2% at the species level and
83.05% at the genus level (taken from Liggenstoffer et al).  For 18S,
settings of 99.0% at the species level and 96.0% at the genus level
seem to produce the most accurate and stable assignments.
