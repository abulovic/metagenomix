from __future__ import absolute_import
import os
from argparse import ArgumentParser

from metagenomix import __version__
from metagenomix.tax import TaxTree
from metagenomix.utils import timeit, get_file_type, get_species_transcript_distribution
from metagenomix.io.sam import parse_cds_sam
from metagenomix.io.megablast import parse_cds_megablast


def get_parse_option_parser():
	parser = ArgumentParser(usage="""parse-mgb <INPUT-FILE> <OUTPUT-DIR> <TYPE>
Parse megablast output file and retrieve semantic alignment data
which is to be stored in the output directory.""",
							version=__version__)
	parser.add_argument(metavar="<INPUT-FILE>", dest="input_file",
		help="Megablast input file (created with -D 3 option)")
	parser.add_argument(metavar="<OUTPUT-DIR>", dest="output_dir",
		help="Output directory to store all the alignment analysis files.")
	parser.add_argument(metavar="<TYPE>", dest="type",
		help="Type of database (cds or genome).")
	return parser

def parse():
	parser = get_parse_option_parser()
	namespace = parser.parse_args()
	if namespace.type != 'cds':
		raise NotImplementedError('Genome parsing will be up in a jiffy.')

	file_type = get_file_type(namespace.input_file)
	with timeit('Parsing alignment file'):
		if file_type == 'blast':
			read_alns, transcripts = parse_cds_megablast(namespace.input_file, namespace.output_dir)
		elif file_type == 'sam':
			read_alns, transcripts = parse_cds_sam(namespace.input_file, namespace.output_dir, binary=False)
		elif file_type == 'bam':
			read_alns, transcripts = parse_cds_sam(namespace.input_file, namespace.output_dir, binary=True)

	with timeit('Loading tax tree'):
		tax_tree = TaxTree()
	spec2trans = get_species_transcript_distribution(transcripts, tax_tree)
	for spec, trans in spec2trans.iteritems():
		if len(filter(lambda t: t.total_coverage > 0.9, trans)) > 0:
			if spec <= 0:
				continue
			int1 = filter(lambda t: t.total_coverage > 0.9, trans)
			print tax_tree.nodes[spec].organism_name,
			print len(trans)
			print ' '.join(map(lambda t: "(%.3f, %.3f)" % (t.total_coverage, t.coverage_fold), int1))
			print
