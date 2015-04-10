from __future__ import absolute_import
import os

import numpy as np
import scipy as sp
from scipy import linalg
from argparse import ArgumentParser

from gensim import corpora

from metagenomix import __version__
from metagenomix.utils import timeit, time_exec
from metagenomix.amino_triplets import nesto, SpeciesCorpus

def get_freq_analysis_parser():
	parser = ArgumentParser(usage="""analyze-prot-freqs <SEQUENCES-DIR> <FREQ-ANALYSIS-TOOL>""",
							version=__version__)
	parser.add_argument(metavar="<SEQ-DIR>", dest="seq_dir",
		help="Directory containing species protein sequences for analysis")
	parser.add_argument(metavar="<FREQ-ANALYSIS-TOOL>", dest="freq_analysis_tool",
		help="Tool for analyzing frequencies in the file")
	return parser

@time_exec('LOADING AA-TRIPS DATA')
def lsi_corpus():
	parser = get_freq_analysis_parser()
	namespace = parser.parse_args()
	seq_dir = namespace.seq_dir
	freq_analysis_tool = namespace.freq_analysis_tool

	with timeit('Calculating species corpus'):
		sc = SpeciesCorpus(seq_dir, namespace.freq_analysis_tool)
		corpora.MmCorpus.serialize('/home/abulovic/tmp/bact_corpus.mm', sc)
		sc.dictionary.save('/home/abulovic/tmp/bact_triplets.dict')

def create_database():
	nesto()