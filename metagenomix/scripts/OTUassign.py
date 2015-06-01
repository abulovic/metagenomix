from __future__ import absolute_import
import os
import time
from argparse import ArgumentParser
from collections import defaultdict
import shutil

from metagenomix import __version__
from metagenomix.tax import TaxTree
import metagenomix.utils as utils
from metagenomix.utils import timeit
from metagenomix.biodb.data_access import DataAccess
import metagenomix.io.blast as blast
import metagenomix.io.sam as sam
import metagenomix.otuassign as otu
import metagenomix.stats.otu as otu_stats
from metagenomix.io.report import Report
import metagenomix.profiling as profiling
from metagenomix.io.config import parse_config
from metagenomix.em import exp_max


def get_OTU_assign_option_parser():
	parser = ArgumentParser(usage="""Args: <INPUT-FILE> <OUTPUT-DIR> <TYPE>
Assign Opearational Taxonomic Units from the BLAST alignment file.
""",
							version=__version__)
	parser.add_argument(metavar="<INPUT-FILE>", dest="input_file",
		help="Blast input file (tabular or xml format)")
	parser.add_argument(metavar="<ORIGINAL-FASTA>", dest="original_fasta",
		help="Original fasta file")
	parser.add_argument(metavar="<OUTPUT-DIR>", dest="output_dir",
		help="Output directory to store all the alignment analysis files.")
	parser.add_argument(metavar="<DB-TYPE>", dest="db_type",
		help="Type of database against which reads were aligned (cds/nt)")
	parser.add_argument('--readcount', dest="read_count", type=int,
		help="Original fasta/fastq read count to for running time estimation")
	return parser

def greedy():
	parser = get_OTU_assign_option_parser()
	args = parser.parse_args()

	if args.read_count != 0:
		total_read_count = args.read_count

	all_reads = utils.reads_from_fasta(args.original_fasta)

	with Report(args.output_dir, 'microbe') as report:

		file_type = utils.get_file_type(args.input_file)

		with report.timeit('Parsing %s <%s>' % (args.input_file, utils.get_appropriate_file_size(args.input_file))):
			if file_type == 'blast':
				count_entries = blast.get_entry_cnt_tab
				parse_func = blast.parse_tab_delimited
			elif file_type in ('sam', 'bam'):
				count_entries = sam.get_entry_cnt_sam
				parse_func = sam.parse_cds_sam
			elif file_type == 'xml':
				count_entries = blast.get_entry_cnt_xml
				parse_func = blast.parse_xml
			else:
				raise ValueError('%s alignment format not supported!' % file_type)

			with utils.timeit('Retrieving entry count'):
				entry_cnt = count_entries(args.input_file)
				report.mark('\tTotal entries: %d' % entry_cnt)
			with utils.timeit('File parsing'):
				read_alns, target_seqs = parse_func(args.input_file, args.db_type, entry_cnt=entry_cnt, detailed=True)
				report.mark('\tTarget sequences: %d' % len(target_seqs))

		with report.timeit('GI2TAX database querying'):
			utils.retrieve_tax_ids(target_seqs, args.db_type)

		with report.timeit('loading tax tree'):
			tt = TaxTree()
			report.mark('Loaded %d nodes.' % len(tt.nodes))

		#exp_max(read_alns, target_seqs, tt)
		#profiling.get_read_overlap(target_seqs, read_alns, tt, all_reads)
		with report.timeit('Generating clusters'):
			clusters = profiling.sequential_read_set_analysis(read_alns, target_seqs, tt)
			report.output_clusters(clusters)

		coverage_limit = 0.6
		fold_limit = 1.
		with report.timeit('greedy transcript assignment'):
			report.mark('Coverage threshold: %.2f' % coverage_limit)
			report.mark('Fold threshold: %.2f' % fold_limit)
			report.mark('#Transcripts (before read assignment)  : %d' % len(target_seqs))
			prefilt_transcripts = otu.greedy_transcript_assign(target_seqs, read_alns)
			report.mark('#Transcripts (after read assignemnt)   : %d' % len(prefilt_transcripts))
			#if args.db_type == 'cds':
			final_transcripts = otu.filter_by_coverage_fold(prefilt_transcripts, 0.0, 1.)
			#else:
			#	final_transcripts = prefilt_transcripts
			report.mark('#Transcripts (after cov-fold filtering): %d' % len(final_transcripts))
			total_reads = len(read_alns)

		with report.timeit('Species assignment stats'):
			s2t_nofilt = utils.get_species_transcript_distribution(target_seqs, tt)
			s2t_greedy = utils.get_species_transcript_distribution(final_transcripts, tt)
			report.mark('#Species (pre-read-assignment) : %d' % len(s2t_nofilt))
			report.mark('#Species (post-read-assignment): %d' % len(s2t_greedy))
			report.rank_distribution(s2t_nofilt, tt, 'nofilt')
			report.rank_distribution(s2t_greedy, tt, 'greedy')
			report.tax_tree(s2t_nofilt, tt, 'nofilt')
			report.tax_tree(s2t_greedy, tt, 'greedy')

		if args.db_type == 'cds':
			new_s2t = otu.remove_orthologue_strains(s2t_greedy)
		else:
			new_s2t = s2t_greedy

		if args.db_type == 'cds':
			with report.timeit('Transcript stats'):
				report.transcript_stats(s2t_nofilt, tt, 'nofilt')
				report.transcript_stats(new_s2t, tt, 'greedy', assigned=True)
				report.tax2reads(new_s2t, tt, 'greedy', 'json')

		if args.db_type == 'cds':
			with report.timeit('Outputing gene expression'):
				report.gene_expression(s2t_nofilt, tt, 'gene_expression_nofilt', assigned=False)
				report.gene_expression(new_s2t, tt, 'gene_expression', assigned=True)

		report.summary(read_alns, s2t_nofilt, new_s2t, args.original_fasta)


def lca():
	parser = get_OTU_assign_option_parser()
	args = parser.parse_args()
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)

	with Report(args.output_dir, 'host-microbe') as report:
		file_type = utils.get_file_type(args.input_file)

		with report.timeit('Parsing %s <%s>' % (args.aln_file, utils.get_appropriate_file_size(args.aln_file))):
			if file_type == 'blast':
				read_alns, target_seqs = blast.parse_tab_delimited(args.input_file, args.db_type)
			elif file_type in ('sam', 'bam'):
				binary = True if file_type == 'bam' else False
				read_alns, target_seqs = sam.parse_cds_sam(args.input_file, args.db_type, binary)
			else:
				raise ValueError('%s alignment format not supported!' % file_type)

		data_access = DataAccess()
		with report.timeit('GI2TAX database querying'):
			gis = set(map(lambda t: t.gi, target_seqs.itervalues()))
			missing_gis = set()
			tax_ids = data_access.get_taxids(gis)
			for target in target_seqs.itervalues():
				target.tax_id = tax_ids.get(target.gi, -1)
				if target.tax_id == -1:
					missing_gis.add(target.gi)
			report.mark('gi count : %d' % len(gis))
			report.mark('tax count: %d' % len(tax_ids))
			report.mark('Number of gis for which no TAXID could be found: %d' % len(missing_gis))

		with report.timeit('Loading tax tree'):
			tt = TaxTree()
			report.mark('Loaded %d nodes.' % len(tt.nodes))

		with report.timeit('lca read assignment'):
			tax2reads = otu.lca_read_assign(target_seqs, read_alns, tt, max_alns=50)
			tax2count = dict(map(lambda i: (i[0], len(i[1])), tax2reads.items()))

		with report.timeit('tax stats'):
			species2reads = utils.get_rank_read_distribution(tax2reads, tt, 'species')
			report.tax2reads(tax2reads, tt, 'lca', 'json')
			report.rank_distribution(species2reads, tt, 'lca')


def get_config_run_parser():
	parser = ArgumentParser(usage="""Args: <metagenomix-CONFIG-FILE>
Run host separation and greedy microbial OTU assignment.
""", version=__version__)
	parser.add_argument(metavar="<metagenomix-CONFIG-FILE>", dest="metagenomix_config",
			help="metagenomix configuration file (in json format)")
	return parser

def host_greedy():
	parser = get_config_run_parser()
	args = parser.parse_args()

	config = parse_config(args.metagenomix_config)
	out_dir = config["output_dir"]

	data_access = DataAccess()

	with Report(out_dir, 'host-microbe') as report:
		if not config["host_separated"]:
			raise utils.NotSupportedError("Host separation currently not supported")

		with report.timeit('loading tax tree'):
			tt = TaxTree()
			report.mark('Loaded %d nodes.' % len(tt.nodes))

		with report.timeit("Loading original fastas"):
			for i, orig_fasta in enumerate(config["input_seq_files"], 1):
				report.load_original_fasta(orig_fasta)

		host_aln_type = config["host_aln_type"].lower()
		for i, aln_file in enumerate(config["host_alignments"]):
			with report.timeit("Extracting read from host alignment"):
				report.extract_reads_from_aln(aln_file, host_aln_type)

		db_type = config["microorganism_db_type"].lower()
		microbe_aln_type = config["microorganism_aln_type"].lower()

		gi_type = 'nucl_gi' if db_type == 'cds' else 'gi'
		get_gi = lambda t: getattr(t, gi_type)

		read_alns = defaultdict(list)
		target_seqs = {}
		for i, aln_file in enumerate(config["microorganism_alignments"]):
			with report.timeit('File parsing'):
				if microbe_aln_type == 'blast':
					read_alns, target_seqs = blast.parse_tab_delimited(aln_file, db_type, 1e5, read_alns, target_seqs,
						filter_low_scoring=False, annotate=False)
				elif microbe_aln_type in ('sam', 'bam'):
					binary = True if microbe_aln_type == 'bam' else False
					read_alns, target_seqs = sam.parse_cds_sam(aln_file, db_type, binary, annotate=False,
						read_alns=read_alns, target_seqs=target_seqs)
				else:
					raise ValueError('%s alignment format not supported!' % microbe_aln_type)

		# Write to file all the rads with alignments to microbes
		report.output_reads(read_alns.keys(), 'microbe', 'reads.txt')
		if microbe_aln_type == 'blast':
			read_alns, target_seqs = blast.annotate_targets(read_alns, target_seqs)
		elif microbe_aln_type in ('sam', 'bam'):
			read_alns, target_seqs = sam.annotate_targets(read_alns, target_seqs)

		with report.timeit('GI2TAX database querying'):
			gis = set(map(lambda t: get_gi(t), target_seqs.itervalues()))
			tax_ids = data_access.get_taxids(gis)
			missing_gis = set()
			for target in target_seqs.itervalues():
				target.tax_id = tax_ids.get(get_gi(target), -1)
				if target.tax_id == -1:
					missing_gis.add(get_gi(target))
			report.mark('gi count : %d' % len(gis))
			report.mark('tax count: %d' % len(set(tax_ids.values())))
			report.mark('Number of gis for which no TAXID could be found: %d' % len(missing_gis))

		with report.timeit('greedy transcript assignment'):
			report.mark('#Transcripts (before read assignment) : %d' % len(target_seqs))
			prefilt_transcripts = otu.greedy_transcript_assign(target_seqs, read_alns)
			report.mark('#Transcripts (after read assignemnt)  : %d' % len(prefilt_transcripts))
			final_transcripts = otu.filter_by_coverage_fold(prefilt_transcripts, 0.5, 1.)
			hypothetical_ids = otu.find_hypothetical(final_transcripts)
			report.mark('#Transcripts (after cov-fold filtering: %d' % len(final_transcripts))
			hypothetical = {_id: final_transcripts[_id] for _id in hypothetical_ids}
			non_hypothetical_ids = set(final_transcripts) - set(hypothetical)
			non_hypothetical = {_id: final_transcripts[_id] for _id in non_hypothetical_ids}

		with report.timeit('Species assignment stats'):
			s2t_nofilt = utils.get_species_transcript_distribution(target_seqs, tt)
			s2t_greedy_nh = utils.get_species_transcript_distribution(non_hypothetical, tt)
			s2t_ss = otu.remove_orthologue_strains(s2t_greedy_nh)
			s2t_greedy_h = utils.get_species_transcript_distribution(hypothetical, tt)
			report.mark('#Species (pre-read-assignment)   : %d' % len(s2t_nofilt))
			report.mark('#Species (post-read-assignment)  : %d' % len(s2t_greedy_nh))
			report.mark('#Species (after-strain-filtering): %d' % len(s2t_ss))
			report.rank_distribution(s2t_nofilt, tt, 'nofilt')
			report.rank_distribution(s2t_ss, tt, 'greedy')
			report.tax_tree(s2t_nofilt, tt, 'nofilt')
			report.tax_tree(s2t_greedy_nh, tt, 'greedy')

		with report.timeit('Transcript stats'):
			report.transcript_stats(s2t_nofilt, tt, 'nofilt')
			report.transcript_stats(s2t_ss, tt, 'greedy', assigned=True)
			report.transcript_stats(s2t_greedy_h, tt, 'hypothetical', assigned=True)

		with report.timeit('Outputing gene expression'):
			report.gene_expression(s2t_nofilt, tt, 'gene_expression/nofilt', assigned=False)
			report.gene_expression(s2t_ss, tt, 'gene_expression/greedy', assigned=True)
			report.gene_expression(s2t_greedy_h, tt, 'gene_expression/hypothetical', assigned=True)

		with report.timeit('Transcript stats'):
			report.transcript_stats(s2t_nofilt, tt, 'nofilt')
			report.transcript_stats(s2t_ss, tt, 'greedy', assigned=True)
			report.tax2reads(s2t_ss, tt, 'greedy', 'json')

		with report.timeit('Different strain stats'):
			report.strains(s2t_ss, tt, 'greedy')
			report.strains(s2t_greedy_h, tt, 'hypothetical')

		with report.timeit('Outputing gene expression'):
			report.gene_expression(s2t_nofilt, tt, 'microbe/gene_expression_nofilt', assigned=False)
			report.gene_expression(s2t_ss, tt, 'microbe/gene_expression', assigned=True)

		total_reads = list()
		for f in os.listdir(out_dir + '/orig_seqs'):
			fname = out_dir + '/orig_seqs/' + f
			with open(fname) as fin:
				total_reads.extend([line.strip() for line in fin])
		total_reads = set(total_reads)
		report.mark('TOTAL')
		report.output_reads(total_reads, 'orig_seqs', 'reads.txt')

		# Let's load all the host-reads
		host_reads = list()
		host_dir = out_dir + os.path.sep + 'host' + os.path.sep
		for f in os.listdir(host_dir):
			with open(host_dir + os.path.sep + f) as fin:
				host_reads.extend([line.strip() for line in fin])
		host_reads = set(host_reads)
		report.mark('HOST')
		report.output_reads(host_reads, 'host', 'reads.txt')

		# Let's load all the microbial reads
		with open(out_dir + '/microbe/reads.txt') as fin:
			microbial_reads = set([line.strip() for line in fin])
		report.mark('MICROBE')
		report.mark('Read count: %d' % len(microbial_reads))

		common_reads = host_reads & microbial_reads
		report.mark('common_reads')
		report.output_reads(common_reads, '', 'common_reads.txt')


		reads_with_no_aln = total_reads - (host_reads | microbial_reads)
		report.mark('NO ALN')
		report.output_reads(reads_with_no_aln, '', 'no_aln_reads.txt')
