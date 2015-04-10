import os
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import seaborn
# import pygal

from metagenomix.stats.otu import get_rank_distribution


def get_read_aln_distribution(read_alns, image_file, original_read_ids=None, original_read_count=None):
	aln_count = np.array(map(lambda alns: len(alns), read_alns.itervalues()))
	_min, _max = aln_count.min(), aln_count.max()

	num_bars = 100 if _max - _min > 100 else _max - _min + 1
	vals, bins, patches = plt.hist(aln_count, bins=_max-_min+1)
	plt.xlabel('Number of alignments')
	plt.ylabel('Number of reads')

	if float(vals[0]) > 0.001 * vals[-1]:
		fig = plt.figure()
		ax1 = fig.add_subplot(211, title='Distribution of alignment count per read')
		ax2 = fig.add_subplot(212)

		max_count = vals[0]
		limit = 0.01 * max_count
		idx = np.argmax(vals < limit)

		ax1.bar(bins[0:idx], vals[0:idx])
		ax1.set_xlabel('Number of alignments')
		ax2.bar(bins[idx+1:], vals[idx:])

		fig.text(0.01, 0.5, 'Number of reads', ha='center', va='center', rotation='vertical')

	plt.tight_layout()
	plt.savefig(image_file)


def get_spec_transcript_distribution(spec2transcript, image_file):
	transcript_count = np.array(map(lambda trans: len(trans), spec2transcript.itervalues()))
	if len(transcript_count) == 0:
		return
	_min, _max = transcript_count.min(), transcript_count.max()
	fig = plt.figure()
	ax = fig.add_subplot(111, title='Number of transcripts per reported species')
	mean = transcript_count.mean()

	num_bars = 100 if _max - _min > 100 else _max - _min + 1
	ax.bar(range(len(transcript_count)), transcript_count)
	ax.hold(True)
	ax.plot([mean]*len(transcript_count), 'r')
	ax.set_xlabel('Species')
	ax.set_ylabel('Number of reported transcripts')
	ax.annotate('mean transcript count = %.2f' % mean, xy=(0, mean), xycoords='data',
				xytext=(0.9, 0.8), textcoords='axes fraction',
                horizontalalignment='right', verticalalignment='top')

	plt.tight_layout()
	plt.savefig(image_file)


def _pie_chart(spec2transcript, tax2color, image, tax_tree, rank='phylum'):
	rank_taxa = defaultdict(int)
	invalid_taxids = (0, -1)
	labels = []
	acceptable_label_percentage = 0.01
	total = 0

	for spec_taxid in spec2transcript.iterkeys():
		tid = tax_tree.get_parent_with_rank(spec_taxid, rank)
		if tid in invalid_taxids:
			continue
		rank_taxa[tid] += 1
		total += 1

	total = float(total)
	counts = []
	for tax, count in sorted(rank_taxa.items(), key=lambda i:i[1]):
		perc = count / total
		if perc >= acceptable_label_percentage:
			labels.append('%s (%.3f)' % (tax_tree.nodes[tax].organism_name, count / total * 100))
		else:
			labels.append('')
		counts.append(count)

	fig = plt.figure()
	ax = fig.add_subplot(111, title='Phylum distribution')
	ax.pie(counts, labels=labels, shadow=True, startangle=90)
	plt.axis('equal')
	plt.savefig(image)
	plt.tight_layout()

def labeled_scatter(spec2transcript, tax_tree, output_dir, naming_prefix, assigned=False):
	# Let's first make a scatter plot of >90% covered transcripts
	def glob(skip):
		i = -1
		l = len(orgs)
		x = 0
		while (True):
			i += 1
			if i <= skip:
				x = yield str(x)
			else:
				if (i - (skip+1)) >= l:
					x = yield 'prst'
				else:
					val = orgs[i - (skip+1)]
					x = yield '\n'.join(val)

	invalid_taxids = (0, -1)
	high_ranking = defaultdict(list)
	for spec, transcripts in spec2transcript.iteritems():
		if spec in invalid_taxids:
			continue
		for trans in transcripts:
			if trans.total_coverage > 0.99 and trans.coverage_fold > 50.:
				high_ranking[spec].append(trans)

	xy_chart = pygal.XY(stroke=False)
	xy_chart.title = 'Transcript coverage-fold'
	xy_chart.x_labels = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	xy_chart.y_labels = [50., 60., 70., 80., 90.]
	skip = len(xy_chart.x_labels) + len(xy_chart.y_labels)
	orgs = []
	for spec, transcripts in high_ranking.iteritems():
		words = tax_tree.nodes[spec].organism_name.split()
		org_name = words[0][0] + '. ' + ' '.join(words[1:])
		values = map(lambda t: (t.total_coverage, t.coverage_fold), transcripts)
		for t in transcripts:
			orgs.append((org_name, t.protein_id))
		xy_chart.add(org_name, values)

	print_fun = glob(skip)
	print_fun.send(None)
	xy_chart.value_formatter = lambda x: print_fun.send(x)
	with open(output_dir + '/pygal_tryout.html', 'w') as fout:
		fout.write(xy_chart.render())

def _scatter_transcript_plot(spec2transcript, tax2color, image,
							 tax_tree, newrank=None, assigned=False):
	colors = []
	coverages = []
	folds = []
	taxids = []
	invalid_taxids = (0, 1)
	for spec_taxid, transcripts in spec2transcript.iteritems():
		if newrank is not None:
			tid = tax_tree.get_parent_with_rank(spec_taxid, newrank)
			if tid in invalid_taxids:
				continue
		else:
			tid = spec_taxid
		for trans in transcripts:
			taxids.append(tid)
			if assigned:
				coverages.append(trans.assigned_total_coverage)
				folds.append(trans.assigned_coverage_fold)
			else:
				coverages.append(trans.total_coverage)
				folds.append(trans.coverage_fold)
			colors.append(tax2color[tid])

	c = Counter(taxids)
	most_common = c.most_common(10)
	mc_colors = defaultdict(list)
	mc_coverages = defaultdict(list)
	mc_folds = defaultdict(list)

	for spec_tax, transcripts in spec2transcript.iteritems():
		tid = spec_tax if newrank is None else tax_tree.get_parent_with_rank(spec_tax, newrank)
		if tid in invalid_taxids:
			continue
		if tid in most_common:
			for trans in transcripts:
				mc_colors[tid].append(tax2color[tid])
				if assigned:
					mc_coverages[tid].append(trans.assigned_total_coverage)
					mc_folds[tid].append(trans.assigned_coverage_fold)
				else:
					mc_coverages[tid].append(trans.total_coverage)
					mc_folds[tid].append(trans.coverage_fold)

	fig = plt.figure()
	ax = plt.subplot(111)

	plt.hold(True)
	for (tid, count) in most_common:
		if tid not in (0, -1):
			plt.scatter(mc_coverages[tid], mc_folds[tid],
						c=tax2color[tid], label=tax_tree.nodes[tid].organism_name)

	plt.scatter(coverages, folds, c=colors)

	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])

	plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
	plt.savefig(image)
	plt.tight_layout()
	plt.hold(False)

def get_transcript_stats(spec2transcript, tax_tree, output_dir, naming_prefix, tax_rank = 'phylum', assigned=False):
	r1 = np.linspace(0, 1, len(spec2transcript))
	species_colors = dict((tax, np.random.rand(3)) for tax in spec2transcript.keys())
	invalid_taxids = (0, -1)

	tax2newrank = {}
	for tax in spec2transcript.iterkeys():
		tax2newrank[tax] = tax_tree.get_parent_with_rank(tax, tax_rank)
	newrank_species = set(tax2newrank.values())
	diff_newrank_species = len(newrank_species)
	r2 = np.linspace(0, 1, diff_newrank_species)
	rank_colors = dict((tax, np.random.rand(3)) for tax in newrank_species)

	spec2count = {k: len(v) for k, v in spec2transcript.iteritems()}
	best10 = list(sorted(spec2count, key=lambda t: spec2count[t], reverse=True))[:10]

	species_scatter_file = '%s%s%s-species.csv' % (output_dir, os.path.sep, naming_prefix)
	with open(species_scatter_file, 'w') as fout:
		fout.write('SpeciesName,Coverage,Fold\n')
		for species, transcripts in spec2transcript.iteritems():
			if species in invalid_taxids:
				org_name = 'unassigned'
			else:
				if species in best10:
					org_name = tax_tree.nodes[species].organism_name
				else:
					org_name = 'Other'

			for trans in transcripts:
				if trans.assigned_total_coverage > 1:
					continue
				fout.write('%s,%.3f,%.3f\n' % (org_name,
											   trans.assigned_total_coverage,
											   trans.assigned_coverage_fold))
	_pie_chart(spec2transcript, rank_colors,
			   '%s%s%s-phylum-piechart.png' % (output_dir, os.path.sep, naming_prefix),
			   tax_tree, 'phylum')

	'''
	_scatter_transcript_plot(spec2transcript, species_colors,
							 os.path.sep.join([output_dir, '%s-species.png' % naming_prefix]), tax_tree,
							 assigned=assigned)
	_scatter_transcript_plot(spec2transcript, rank_colors,
							 os.path.sep.join([output_dir, '%s-%s.png' % (naming_prefix, tax_rank)]), tax_tree, tax_rank,
							 assigned=assigned)
	'''
	#labeled_scatter(spec2transcript, tax_tree, output_dir, naming_prefix, assigned)
