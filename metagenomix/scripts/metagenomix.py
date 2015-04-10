#! /usr/bin/env python
# encoding: utf-8

"""metagenomix [options] <FASTQ FILE> <OUTPUT_DIR>

Metagenomic analysis pipeline is divided in steps.
Each of the steps can be skipped it the intermediate result is already present.
1. Preprocessing
2. Assembly (optional, default=False)
3. Alignment
4. Binning
The intermediate steps are stored in the output directory.
Since any and all of these tasks can be extremely time-consuming,
if you already have a generated intermediate result for any of the
steps, you can provide it via command-line options.
The metagenomic binning analysis can be performed in two ways:
1) Simple mode, where only taxonomic assignments of the
   alignments are considered.
2) Elaborate mode, where the 16S content is analyzed and the
   the genomic region coverate is used to estimate the
   sample evenness and complexity.

"""
from __future__ import absolute_import

import os
from argparse import ArgumentParser

from metagenomix import __version__
from metagenomix.preprocessing import run_preprocessing
from metagenomix.graph import PipelineExecGraph


def get_option_parser():
	parser = ArgumentParser(usage=__doc__, version=__version__)

	parser.add_argument(metavar="<FASTQ_FILE>", dest="fastq_file",
		help="Metagenomic / Metatranscriptomic fastq file")
	parser.add_argument(metavar="<OUTPUT_DIR>", dest="output_dir",
		help="Output directory used for storing intermediate and final results.")
	parser.add_argument("--single-stage", action="store", metavar="SINGLE_STAGE", dest="single_stage",
		choices=['preprocessing', 'alignment', 'binning'],
		help="Run only a single stage of the metagenomic pipeline. Options: (preprocessing, alignment, binning)")
	parser.add_argument("--override", action="store_true", dest="override", default=False,
		help="If supplied, the pipeline will override the steps already performed and perform them again.")

	group = parser.add_argument_group("Preprocessing options",
		description="Options which tune the preprocessing step of the metagenomix.")
	group.add_argument("--fastqc", action="store", metavar="FASTQC_REPORT", dest="fastqc_report",
		help="Path to the already generated FastQC report.")
	group.add_argument('--contaminants', action="store", metavar="CONTAMINANTS", dest="contaminants",
		help="Path to the contaminants file (used by FastQC)")

	group = parser.add_argument_group("Alignment options",
		description="Options to tune the alignment step of the metagenomix.")
	group.add_argument("--blast_db", action="store", metavar="BLAST_DB", dest="blast_db",
		help="Path to the blast database to use for the alignment step.")
	group.add_argument("--aln_file", action="store", metavar="ALN_FILE", dest="aln_file",
		help="Already generated blast alignment file (default: megablast -D 3 output option)")
	group.add_argument("--blast_format", action="store", metavar="BLAST_FORMAT", dest="blast_format",
		help="Blast format used for parsing alignment file if one is provided, or for generating if not.")

	return parser


def main(cmdlineargs=None):

	parser = get_option_parser()
	namespace = parser.parse_args()
	assert(os.path.isdir(namespace.output_dir))
	if namespace.single_stage is not None:
		namespace.single_stage = namespace.single_stage.capitalize()
	execution_graph = PipelineExecGraph(namespace.output_dir, stage=namespace.single_stage, override=namespace.override)

	execution_graph.nodes['Preprocessing'].task = run_preprocessing
	execution_graph.execute(args=namespace)

if __name__ == '__main__':
	main()
