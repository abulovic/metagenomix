from __future__ import absolute_import
import os
import time
from argparse import ArgumentParser
from collections import defaultdict
from contextlib import nested

from Bio import SeqIO

from metagenomix import __version__
from metagenomix.tax import TaxTree
import metagenomix.utils as utils
from metagenomix.biodb.data_access import DataAccess
import metagenomix.io.db as db_utils
import metagenomix.io.out as out

def get_extract_reads_parser():
	parser = ArgumentParser(usage="""Extract reads from a fasta file """,
							version=__version__)
	parser.add_argument(metavar="<READS>", dest="reads",
		help="File with newline separated list of read identifiers")
	parser.add_argument(metavar="<ORIGINAL-FASTA>", dest="orig_fasta",
		help="Original fasta")
	parser.add_argument(metavar="<OUTPUT-FILE>", dest="output_file",
		help="Output fasta file with reads specified in the reads file")
	return parser

def get_extract_subtaxa_parser():
	parser = ArgumentParser(usage="""Extract all taxa from a database that are subtaxa of the specified parent tax.""",
							version=__version__)
	parser.add_argument(metavar="<INPUT-DB>", dest="input_db",
		help="Input database (format as downloaded from NCBI site)")
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('--tax', type=int, dest="parent_tax",
		help="Parent tax for which all children taxa will be included")
	group.add_argument('--file', type=str,
		help="Text file with a single line of comma-separated parent taxa")
	parser.add_argument(metavar="<OUTPUT-DIR>", dest="output_dir",
		help="Path to output database file")
	parser.add_argument('--merge', action="store_true",
		help="Merge all species records in a single file")
	parser.add_argument('--gene-count', type=int, dest='gene_count',
		help="Number of non-hypothetical different genes to choose from each species")
	parser.add_argument('--db-type', type=str, dest="db_type", default="cds",
		help="Database type for header parsing", choices=("cds", "nt", "genome"))
	return parser

def extract_reads():
	parser = get_extract_reads_parser()
	args = parser.parse_args()
	file_type = utils.get_file_type(args.orig_fasta)

	with open(args.reads) as fin:
		reads = [l.strip() for l in fin]
	out.extract_reads(reads, args.orig_fasta, args.output_file, file_type)


def extract_subtaxa():
	parser = get_extract_subtaxa_parser()
	args = parser.parse_args()

	if args.db_type == "cds":
		parse_header = db_utils.parse_cds_header
		def get_taxid(data, data_access):
			return int(data["taxon"])
	elif args.db_type == "genome":
		parse_header = db_utils.parse_genome_header
		def get_taxid(data, data_access):
			taxid = data_access.get_taxids((int(data["gi"]),), format=list)
			if len(taxid) == 0:
				return None
			else:
				return taxid[0]
	else:
		parse_header = db_utils.parse_nt_header

	if args.parent_tax:
		taxa = [args.parent_tax]
	elif args.file:
		with open(args.file) as fin:
			taxa = map(lambda t: int(t), fin.read().strip().split(','))

	tt = TaxTree()

	if args.merge:
		single_file = '{0}{1}db-extract.fa'.format(args.output_dir, os.path.sep)
		handle = open(single_file, 'w')
		write = lambda tax, record: SeqIO.write(record, handle, 'fasta')
		close = lambda: handle.close()
	else:
		tax2handle = {}
		for tax in taxa:
			fname = '{0}{1}{2}.fa'.format(args.output_dir, os.path.sep, tt.get_org_name(tax))
			tax2handle[tax] = open(fname, 'w')
			write = lambda tax, record: SeqIO.write(record, tax2handle[tax], 'fasta')
			close = lambda: map(lambda h: h.close(), tax2handle.values())

	data_access = DataAccess()


	if args.gene_count is None:
		with open(args.input_db, 'r') as fin:
			records = SeqIO.parse(fin, 'fasta')
			for record in records:
				data = parse_header(record.name)
				taxid = get_taxid(data, data_access)
				if taxid is None:
					continue
				for pt in taxa:
					if tt.is_child(taxid, pt) or taxid == pt:
						write(pt, record)
	else:
		tax2genes = defaultdict(set)
		interesting_taxa = set(taxa)
		with open(args.input_db) as fin:
			records = SeqIO.parse(fin, 'fasta')
			for record in records:
				if not interesting_taxa:
					break
				data = parse_header(record.name)
				taxid = get_taxid(data, data_access)
				if taxid is None:
					continue
				lineage = tt.get_lineage(taxid)
				overlap = set(lineage) & interesting_taxa
				if not overlap:
					continue
				parent_tax = overlap.pop()
				if 'product' not in data:
					continue
				product = data['product']
				if product != 'hypothetical protein' and product not in tax2genes[parent_tax]:
					write(parent_tax, record)
					tax2genes[parent_tax].add(product)
					if len(tax2genes[parent_tax]) >= args.gene_count:
						interesting_taxa.remove(parent_tax)
	close()
