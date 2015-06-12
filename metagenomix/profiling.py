from contextlib import nested
from collections import defaultdict

import numpy as np

import metagenomix
import metagenomix.io.db as db_utils
import metagenomix.utils as utils
from metagenomix.io.blast import parse_tab_delimited, parse_xml
from metagenomix.io.sam import parse_cds_sam
from metagenomix.biodb.data_access import DataAccess


class OrgCluster(object):
	def __init__(self, lca, strains, tax_tree, ureads, ireads, total_reads):
		self.lca = lca
		self.strains = strains
		self.tax_tree =  tax_tree
		self.ureads = ureads
		self.ireads = ireads
		self.total_reads = total_reads
		self.percentage = float(len(self.ureads)) / self.total_reads

	def get_name(self):
		return self.tax_tree.get_org_name(self.lca)

	def get_rank(self):
		return self.tax_tree.get_org_rank(self.lca)

	def get_percentage(self):
		return self.percentage

	def __str__(self):
		return '%s|%s|%.2f|%d' % (self.get_name(),
							      self.get_rank(),
							      self.get_percentage()*100,
							      len(self.ureads))

	def __lt__(self, other):
		return len(self.ureads) < len(other.ureads)


def tax_to_reads(read_alns, target_seqs, tax_tree, filter_low_cnt=False):
	species2reads = defaultdict(set)
	for read, alns in read_alns.iteritems():
		for aln in alns:
			tid = target_seqs[aln.target_id].tax_id
			species2reads[tid].add(read)
	if filter_low_cnt:
		low_cnt_taxa = set()
		for tax, reads in species2reads.iteritems():
			if len(reads) <= 5:
				low_cnt_taxa.add(tax)
		for tax in low_cnt_taxa:
			del species2reads[tax]

	return species2reads


def get_simple_profile(alignment_file, aln_type, db_type, tax_tree):
	assert(aln_type in metagenomix.supported_aln_types)
	assert(db_type in metagenomix.supported_db_types)

	if aln_type == 'blast':
		read_alns, target_seqs = parse_tab_delimited(alignment_file, db_type, annotate=False)
	elif aln_type == 'sam':
		read_alns, target_seqs = parse_cds_sam(alignment_file, db_type, annotate=False)
	elif aln_type == 'xml':
		read_alns, target_seqs = parse_xml(alignment_file, db_type, annotate=False)
	utils.retrieve_tax_ids(target_seqs, db_type)

	return tax_to_reads(read_alns, target_seqs, True)

def get_read_overlap(target_seqs, read_alns, tax_tree, all_reads):
	tax2reads = tax_to_reads(read_alns, target_seqs, tax_tree, True)
	tax2target = {ts.tax_id: ts for ts in target_seqs.itervalues()}
	target_cnt = len(tax2reads)
	intersect_matrix = np.zeros((target_cnt, target_cnt), np.int32)
	diff_matrix = np.zeros((target_cnt, target_cnt), np.int32)
	sorted_taxa = sorted(tax2reads, reverse=True, key=lambda k: len(tax2reads[k]))
	targets = [tax2target[tax] for tax in sorted_taxa]
	print len(targets)
	limit = min(50, len(sorted_taxa))
	for i, tax1 in enumerate(sorted_taxa[:limit]):
		for j, tax2 in enumerate(sorted_taxa[:limit]):
			intersect_matrix[i,j] = len(tax2reads[tax1] & tax2reads[tax2])
			diff_matrix[i,j] = len((all_reads - tax2reads[tax1]) & tax2reads[tax2])
	with nested(open('intersect_matrix.txt', 'w'), open('diff_matrix.txt', 'w')) as (ifout, dfout):
		print intersect_matrix.shape
		#for i in xrange(intersect_matrix.shape[0]):
		for i in range(limit):
			ifout.write('%s (%d, %.4f)\t' % (tax_tree.get_org_name(sorted_taxa[i]), targets[i].alignment.length, targets[i].alignment.get_coverage()))
			ifout.write(' '.join([str(entry) for entry in intersect_matrix[i]][:limit]))
			ifout.write('\n')
			dfout.write('%s (%d, %.4f)\t' % (tax_tree.get_org_name(sorted_taxa[i]), targets[i].alignment.length, targets[i].alignment.get_coverage()))
			dfout.write(' '.join([str(entry) for entry in diff_matrix[i]][:limit]))
			dfout.write('\n')

def choose_strain(tax2target, cluster):
	targets = []
	target_probs = []
	total_len = 0
	total_fold = 0.
	for tax in tax2target.keys():
		if tax in cluster:
			target = tax2target[tax]
			targets.append(target)
			total_len += target.alignment.length
			total_fold += target.alignment.get_fold()
	for target in targets:
		#target_probs.append(target.alignment.length * target.alignment.get_coverage()\
		#					* target.alignment.get_fold() / (total_len * total_fold))
		target_probs.append(target.alignment.get_coverage())

	target2prob = {t:p for t, p in zip(targets, target_probs)}
	return sorted(target2prob.keys(), reverse=True, key=lambda t: target2prob[t])[0]

def sequential_read_set_analysis(read_alns, target_seqs, tax_tree):
	aln_reads = set(read_alns.keys())
	tax2reads = {ts.tax_id: ts.alignment.reads for ts in target_seqs.values()}
	#tax2reads = tax_to_reads(read_alns, target_seqs, tax_tree, True)
	tax2target = {ts.tax_id: ts for ts in target_seqs.itervalues()}
	sorted_taxa = sorted(tax2reads, reverse=True, key=lambda k: len(tax2reads[k]))

	look_further = True
	cluster_reads = set()
	cluster_taxa = set()
	new_cluster = False
	total_used_reads = set()
	subset_taxa = set()
	tclusters = []
	runion = []
	rintersect = []

	for tax in sorted_taxa:
		if len(total_used_reads) >= 0.999*len(aln_reads):
			print 'FINISHED'
			break
		reads = tax2reads[tax]
		if not tclusters:
			tclusters.append(set([tax]))
			runion.append(set(reads))
			rintersect.append(set(reads))
			total_used_reads |= reads
		else:
			solved = False
			for i, (tc, rc) in enumerate(zip(tclusters, runion)):
				#spec_tax = []
				#for cluster_tax in tc:
				#	spec_tax.append(tax_tree.get_parent_with_rank(cluster_tax, 'species'))
				#current_spec_tax = tax_tree.get_parent_with_rank(tax, 'species')
				#if current_spec_tax in spec_tax:
				#	tc.add(tax)
				#	rc |= reads
				#	total_used_reads |= reads
				#	continue
				overlap = reads & rc
				if len(overlap) > 0.95 * len(rc):
					tc.add(tax)
					rc |= reads
					rintersect[i] &= reads
					total_used_reads |= reads
					solved = True
				else:
					if len(overlap) >= 0.95 * len(reads):
						subset_taxa.add(tax)
						solved = True

			if not solved:
				unique = reads - total_used_reads
				if len(unique) == 0:
					subset_taxa.add(tax)
				else:
					tclusters.append({tax})
					runion.append(set(unique))
					rintersect.append(set(unique))
					total_used_reads |= unique

	# postprocessing:
	print 'Clusters before postprocessing:', len(tclusters)
	i_to_delete = []
	for i, tc in enumerate(tclusters):
		len_intersect = len(rintersect[i])
		len_union = len(runion[i])
		if len_intersect < 0.7 * len_union:
			i_to_delete.append(i)
	i_to_delete.sort(reverse=True)
	for i in i_to_delete:
		tclusters.pop(i)
		rintersect.pop(i)
		runion.pop(i)
	print 'Clusters after postprocessing :', len(tclusters)
	print

	total_perc = 0
	total_reads = len(read_alns)
	clusters = []
	for i, r in sorted(enumerate(runion), reverse=True, key=lambda v: len(v[1])):
		if len(runion[i]) == 0:
			continue
		lca = tax_tree.find_lca(tclusters[i])
		clusters.append(OrgCluster(lca, tclusters[i], tax_tree, r, rintersect[i], total_reads))
		print 'Cluster[%2d]: %s (%s)' % (i, tax_tree.get_org_name(lca), tax_tree.get_org_rank(lca))
		intersect_len = len(rintersect[i])
		union_len = len(runion[i])
		print '#READS (UNION):    %d/%d' % (union_len, total_reads)
		print '#READS (INTSC):    %d/%d' % (intersect_len, total_reads)
		for tax_id in tclusters[i]:
			print tax_tree.get_org_name(tax_id), '(%.2f)' % (tax2target[tax_id].alignment.get_coverage()*100), ' (%d)' % tax2target[tax_id].alignment.length
		total_perc += union_len / float(total_reads)
		print 'Cumsum (perc):  ', total_perc
		print

	for c in clusters:
		c.percentage = c.get_percentage() / total_perc

	# Let's remove the duplicates now
	tax2clusters = defaultdict(list)
	for c in clusters:
		tax2clusters[c.lca].append(c)
	new_clusters = []
	for tax, clusters in tax2clusters.iteritems():
		if len(clusters) == 1:
			new_clusters.append(clusters[0])
		else:
			lca = clusters[0].lca
			ureads = set()
			ireads = set()
			strains = set()
			for c in clusters:
				ureads |= c.ureads
				ireads &= c.ireads
				strains |= c.strains
			new_clusters.append(OrgCluster(lca, strains, tax_tree, ureads, ireads, total_reads))

	new_clusters.sort(reverse=True)

	# Remove all clusters with 0 unique reads
	i_to_delete = []
	for i, c in enumerate(new_clusters):
		unique = 0
		for read_id in c.ureads:
			if len(read_alns[read_id]) == 1:
				unique += 1
		total = len(c.ureads)
		if tax_tree.get_org_rank(c.lca) in ('class', 'phylum', 'kingdom', 'superkingdom') or\
		   c.lca in (1,131567):
			print 'BAD CLUSTER', i, c, tax_tree.get_org_rank(c.lca)
			max_overlap = 0
			max_overlap_clust = None
			for clust in new_clusters:
				if clust.lca == c.lca:
					continue
				ureads = reduce(lambda a, b: a|b, [tax2reads[s] for s in clust.strains])
				overlap = len(c.ureads & ureads)
				if overlap > max_overlap:
					max_overlap = overlap
					max_overlap_clust = clust
			i_to_delete.append(i)
			if max_overlap_clust is not None:
				max_overlap_clust.ureads |= c.ureads

	print 'to delete:', i_to_delete
	for i in reversed(i_to_delete):
		new_clusters.pop(i)

	# Remove all clusters with 0 unique reads
	i_to_delete = []
	for i, c in enumerate(new_clusters):
		unique = 0
		for read_id in c.ureads:
			if len(read_alns[read_id]) == 1:
				unique += 1
		total = len(c.ureads)
		#import pdb
		#pdb.set_trace()
		if total <= 5 or unique == 0 or (unique/float(total) < 0.05 and unique < 5):
			max_overlap = 0
			max_overlap_clust = None
			for clust in new_clusters:
				if clust.lca == c.lca:
					continue
				ureads = reduce(lambda a, b: a|b, [tax2reads[s] for s in clust.strains])
				overlap = len(c.ureads & ureads)
				if overlap > max_overlap:
					max_overlap = overlap
					max_overlap_clust = clust
			if max_overlap_clust is None:
				if total <= 5:
					i_to_delete.append(i)
				continue
			i_to_delete.append(i)
			max_overlap_clust.ureads |= c.ureads

	for i in reversed(i_to_delete):
		new_clusters.pop(i)

	# TODO! MAKE THIS BETTER!
	union_len = sum([len(c.ureads) for c in new_clusters])
	total_perc = union_len / float(total_reads)
	for c in new_clusters:
		c.percentage = len(c.ureads) / float(union_len)

	return new_clusters
