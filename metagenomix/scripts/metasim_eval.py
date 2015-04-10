from __future__ import absolute_import
import os
from argparse import ArgumentParser

from metagenomix import __version__
from metagenomix.testing.metasim import get_read_count, synth_dataset_profile
from metagenomix.tax import TaxTree
import metagenomix.io.out as out

def get_msim_profile_parser():
	parser = ArgumentParser(usage="""Args: <METASIM-FASTA> <OUTPUT-JSON>
Create a species -> targets -> reads json from metasim fasta file.""")
	parser.add_argument(metavar='<METASIM-FASTA>', dest='metasim_fasta')
	parser.add_argument(metavar='<OUTPUT-JSON>', dest='output_json')
	return parser

def get_metasim_eval_parser():
	parser = ArgumentParser(usage="""Args: <metagenomix-JSON> <METASIM-FASTA> <OUTPUT-DIR>
Create an metagenomix performance evaluation from MetaSim data.
""",
							version=__version__)
	parser.add_argument(metavar="<metagenomix-SPECIES-JSON>", dest="metagenomix_json",
		help="metagenomix generated JSON with species tax: read ID list")
	parser.add_argument(metavar="<METASIM-FASTA>", dest="metasim_fasta",
		help="Fasta file created using MetaSim.")
	parser.add_argument(metavar="<OUTPUT-DIR>", dest="output_dir",
		help="A directory where the analysis data will be stored")
	return parser

def metasim_analysis():
	parser = get_metasim_eval_parser()
	args = parser.parse_args()

	tt = TaxTree()
	metasim_json = '{0}{1}metasim-species-reads.json'.format(args.output_dir, os.path.sep)
	synth_dataset_profile(args.metasim_fasta, tt, metasim_json)

	reads_comp_json = '{0}{1}metagenomix-metasim-reads-comparison.json'.format(args.output_dir, os.path.sep)
	out.read_comparison_json(args.metagenomix_json, metasim_json, reads_comp_json, tt)
	transcript_comp_json = '{0}{1}metagenomix-metasim-transcript-comparison.json'.format(args.output_dir, os.path.sep)
	out.transcript_comparison_json(args.metagenomix_json, metasim_json, transcript_comp_json, tt)

	template_html = '/home/abulovic/BINNER/tools/metagenomix/metagenomix/visualization/metasim-metagenomix-template.html'
	reads_output_html = '{0}{1}metagenomix-metasim-reads-comparison.html'.format(args.output_dir, os.path.sep)
	out.default_d3_plot(reads_comp_json, template_html, reads_output_html, 2)

	transcripts_output_html = '{0}{1}metagenomix-metasim-transcript-comparison.html'.format(args.output_dir, os.path.sep)
	out.default_d3_plot(transcript_comp_json, template_html, transcripts_output_html, 2)


def msim_profile():
	parser = get_msim_profile_parser()
	args = parser.parse_args()

	tt = TaxTree()
	synth_dataset_profile(args.metasim_fasta, tt, args.output_json)
