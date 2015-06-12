from __future__ import absolute_import
import os
import time
import shutil
from collections import defaultdict
from contextlib import contextmanager

from Bio import SeqIO

from metagenomix import __version__
import metagenomix.stats.aln as aln_stats
import metagenomix.stats.otu as otu_stats
import metagenomix.utils as utils
import metagenomix.io.blast as blast
import metagenomix.io.sam as sam
import metagenomix.io.out as out
import metagenomix.visualization as viz


class Report(object):

	def __init__(self, output_dir, dir_type='host-microbe'):
		super(Report, self).__init__()
		self.output_dir = output_dir
		self.create_dir_structure(dir_type)
		fname = '%s%sREPORT.txt' % (self.output_dir, os.path.sep)
		self.report_file = open(fname, 'w')

	def __enter__(self):
		self._report('metagenomix version %s' % __version__)
		date = time.strftime("%d/%m/%Y")
		_time = time.strftime("%H:%M:%S")
		self._report('Date: %s, time: %s' % (date, _time), 2)
		self._report('Output folder: %s' % self.output_dir, 2)
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		self.report_file.close()

	def _report(self, message, ns=1):
		print message, '\n' * (ns - 1)
		self.report_file.write('%s%s' % (message, '\n'*ns))

	def mark(self, message, ns=1):
		self._report(message, ns)

	@contextmanager
	def timeit(self, message):
		self._report('Started <%s>' % message)
		start = time.time()
		yield
		stop = time.time()
		duration = stop - start
		self._report('<%s> ended: execution lasted %.3f sec.' % (message, duration), 2)

	def create_dir_structure(self, dir_type='host-microbe', erase=True):
		if dir_type not in ('host-microbe', 'microbe'):
			raise utils.NotSupportedError('Cannot create output directory structure of type "%s"', dir_type)
		utils.create_dir(self.output_dir, erase)
		if dir_type == 'host-microbe':
			utils.create_dir('%s/host' % self.output_dir, erase)
			utils.create_dir('%s/orig_seqs' % self.output_dir, erase)
			utils.create_dir('%s/microbe' % self.output_dir, erase)
			utils.create_dir('%s/plots' % self.output_dir, erase)
			utils.create_dir('%s/gene_expression' % self.output_dir, erase)
			utils.create_dir('%s/composition' % self.output_dir, erase)
		elif dir_type == 'microbe':
			utils.create_dir('%s/composition' % self.output_dir, erase)
			utils.create_dir('%s/gene_expression' % self.output_dir, erase)
			utils.create_dir('%s/plots' % self.output_dir, erase)
		c3_css = os.path.join(viz.__path__[0], 'c3.css')
		shutil.copyfile(c3_css, '%s/plots/c3.css' % (self.output_dir))


	def load_original_fasta(self, fasta):
		fsize = utils.get_appropriate_file_size(fasta)
		self._report("%s <%s>" % (fasta, fsize))
		reads = set()
		with open(fasta) as fin:
			records = SeqIO.parse(fin, 'fasta')
			for rec in records:
				reads.add(rec.name.split()[0])
		self._report("  Read count: %d" % len(reads))
		fname = '.'.join(fasta.split(os.path.sep)[-1].split('.')[:-1]) + '.txt'
		out.output_reads(reads, self.output_dir + os.path.sep + 'orig_seqs', fname, sep='\n')

	def output_reads(self, reads, subdir, fname):
		self._report("Read count: %d" % len(reads))
		out_dir = os.path.sep.join([self.output_dir, subdir])
		out.output_reads(reads, out_dir, fname)

	def output_clusters(self, clusters):
		out_file = os.path.join(self.output_dir, 'composition', 'clusters.txt')
		out.output_clusters(clusters, out_file)

	def extract_reads_from_aln(self, alignment_file, aln_type, subdir='host'):
		fsize = utils.get_appropriate_file_size(alignment_file)
		self._report("ALN FILE: %s <%s>" % (alignment_file, fsize))
		if aln_type == 'blast':
			reads = blast.extract_unique_reads(alignment_file)
		elif aln_type in ('sam', 'bam'):
			reads = sam.extract_unique_reads(alignment_file)
		self._report("  Unique read count: %d" % len(reads))
		out_dir = os.path.sep.join([self.output_dir, subdir])
		fname = '.'.join(alignment_file.split(os.path.sep)[-1].split('.')[:-1]) + '.txt'
		out.output_reads(reads, out_dir, fname)

	def transcript_stats(self, species2transcript, tax_tree,
							   naming_prefix, assigned=False,
							   additional_ranks=('phylum',)):
		# TODO: correct the aln_stats method to receive multiple ranks
		aln_stats.get_transcript_stats(species2transcript, tax_tree, self.output_dir,
									   naming_prefix, tax_rank='phylum', assigned=assigned)
		json_file = '{0}{1}composition{1}species-transcript-reads.json'.format(self.output_dir, os.path.sep)
		out.spec_transcript_read_json(species2transcript, tax_tree, json_file)

	def gene_expression(self, species2transcript, tax_tree, prefix, assigned=False):
		out_dir = self.output_dir + os.path.sep + prefix
		out.output_species_transcripts(species2transcript, out_dir, tax_tree, assigned=assigned)
		files = out.json_transcripts(species2transcript, out_dir, tax_tree)
		template_html = os.path.join(viz.__path__[0], 'scatter-template.html')
		for (json_covfilt, json_names) in files:
			name = json_covfilt.split(os.path.sep)[-1].split('.')[0].split('-')[0]
			output_html = '{0}{1}plots{1}{2}-scatter.html'.format(self.output_dir,
					       os.path.sep, name)
			out.c3_scatter(json_covfilt, json_names, template_html, output_html)

	def rank_distribution(self, species2transcript, tax_tree,
								naming_prefix, additional_ranks=('phylum',)):
		out_dir = '%s%scomposition' % (self.output_dir, os.path.sep)
		_, files = otu_stats.simple_taxa_transcript_distribution(species2transcript, out_dir,
													  		  tax_tree, naming_prefix, additional_ranks)
		template = os.path.join(viz.__path__[0], 'piechart-template.html')
		for f in files:
			name = f.split(os.path.sep)[-1].split('.')[-2]
			output_html = '{0}{1}plots{1}{2}-piechart.html'.format(self.output_dir,
					       os.path.sep, name)
			out.default_d3_plot(f, template, output_html)
		# finaly, let's output species to read count info
		fname = '{0}{1}composition{1}{2}-readcount.csv'.format(self.output_dir,
															   os.path.sep,
															   naming_prefix)
		#out.species_readcount(species2transcript, tax_tree, fname)


	def tax2reads(self, tax2reads, tax_tree, naming_prefix, out_format):
		fname = '%s%s%s-tax2reads.%s' % (self.output_dir, os.path.sep,
										 naming_prefix, out_format)
		self._report('Writing tax2read stats to:\n%s' % fname)

	def tax_tree(self, species2transcript, tax_tree, naming_prefix):
		fname = '{0}{1}composition{1}{2}-tree.json'.format(self.output_dir, os.path.sep,
											 naming_prefix)
		template = os.path.join(viz.__path__[0], 'dendrogram-template.html')
		output_html = '{0}{1}plots{1}{2}tax-tree.html'.format(self.output_dir, os.path.sep, naming_prefix)
		if naming_prefix == 'clusters':
			out.tax_tree_from_clusters(species2transcript, tax_tree, fname)
		else:
			out.output_tax_tree(species2transcript, tax_tree, fname)
		out.default_d3_plot(fname, template, output_html)

	def strains(self, species2transcript, tax_tree, prefix):
		fname = os.path.sep.join([self.output_dir, 'microbe', '%s-strains.txt' % prefix])
		with open(fname, 'w') as fout:
			for (spec_tax, transcripts) in sorted(species2transcript.items(), reverse=True, key=lambda i: len(i[1])):
				if spec_tax in (0, -1):
					continue
				strain2count = defaultdict(int)
				fout.write('%s\n' % tax_tree.nodes[spec_tax].organism_name)
				for t in transcripts:
					strain2count[t.tax_id] += 1
				for strain_tax, count in sorted(strain2count.items(), reverse=True, key=lambda i: i[1]):
					fout.write('%s, %d\n' % (tax_tree.nodes[strain_tax].organism_name, count))

	def summary(self, read_aln, species_nofilt, species_final, input_fasta):
		summary_data = {}
		read_cnt = utils.get_seq_count(input_fasta)
		summary_data['Reads in original sample'] = read_cnt

		reads_with_aln_cnt = len(read_aln)
		summary_data['Reads with at least one alignment'] = reads_with_aln_cnt

		avg_aln_per_reads = sum([len(alns) for alns in read_aln.itervalues()]) / float(reads_with_aln_cnt)
		summary_data['Average number of alignments per read'] = avg_aln_per_reads

		all_species = len(species_nofilt)
		summary_data['Species with at least one reported alignment'] = all_species

		reported_species = len(species_final)
		summary_data['Final reported species'] = reported_species

		summary_file = os.path.join(self.output_dir, 'run-summary.txt')
		with open(summary_file, 'w') as fout:
			for name, value in summary_data.iteritems():
				fout.write('{}:{}\n'.format(name, value))
