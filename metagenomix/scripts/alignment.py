from __future__ import absolute_import

import os
import sqlite3
from argparse import ArgumentParser

from metagenomix import __version__
from metagenomix.tax import TaxTree
from metagenomix.utils import table_exists, timeit
import metagenomix.aln_analysis as analysis
from metagenomix.io.datainput import save_reads_to_fasta, export_OTU_to_csv, export_OTU_to_json

def get_aln_option_parser():
	parser = ArgumentParser(usage="""analyze-alignment <DB-FILE>
Utility to analyze the alignment from the database.
Statistics provided are number of reads with no alignments,
histogram of number of alignments per read etc.""",
							version=__version__)
	parser.add_argument(metavar="<DB-FILE>", dest="db_file",
		help="Database file where the read table will be stored")
	parser.add_argument(metavar="<OUTPUT-DIR>", dest="output_dir",
		help="Output directory to store all the analysis files.")
	return parser


def analyze():
	parser = get_aln_option_parser()
	namespace = parser.parse_args()
	db = sqlite3.connect(namespace.db_file)
	cursor = db.cursor()
	empty_reads_file = os.path.sep.join([namespace.output_dir, 'no_alignment_reads.fa'])
	basic_stats_file = os.path.sep.join([namespace.output_dir, 'basic_alignment_stats.txt'])
	phylum_composition_file = os.path.sep.join([namespace.output_dir, 'phylum_composition.csv'])
	no_otu_reads_file = os.path.sep.join([namespace.output_dir, 'no_OTU_assignment_reads.fa'])
	otu_csv = os.path.sep.join([namespace.output_dir, 'OTU_assignment.csv'])
	otu_json = os.path.sep.join([namespace.output_dir, 'OTU_assignment.json'])

	with timeit('LOADING TAX TREE'):
		tt = TaxTree()
	read2tax, low_count_taxa = analysis.perform_OTU_analysis_on_db(db, cursor, tt, empty_reads_file, basic_stats_file, phylum_composition_file)
	taxa = set(read2tax.values())
	num_alns = len(read2tax)

	no_otu_reads = filter(lambda r: read2tax[r] == -1, read2tax.keys())
	save_reads_to_fasta(db, cursor, no_otu_reads, no_otu_reads_file)

	export_OTU_to_csv(read2tax, tt, otu_csv)
	export_OTU_to_json(read2tax, tt, otu_json)

	db.close()
