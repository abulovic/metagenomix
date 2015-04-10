'''
A module for assigning Operational Taxonomic Unit from a set
of transcripts with information on coverage, number of reads
aligned and the like.
'''
from collections import defaultdict, deque

from metagenomix.utils import simple_cds_orthologue_test, get_species_transcript_distribution



def assign_OTU_from_transcripts(transcripts, tax_tree):
	spec2transcript = get_species_transcript_distribution(transcripts, tax_tree)

def greedy_transcript_assign(transcripts, read_aln, filter_low_coverage=False):
	l1 = len(transcripts)
	if filter_low_coverage:
	# TODO: check how many transcripts we just loose using this step
		transcripts = dict(filter(lambda t: len(transcripts[t[0]].get_alignments()) > 5, transcripts.iteritems()))
		l2 = len(transcripts)
	# We first sort transcripts by their coverage and fold
	assigned = set()
	sorted_transcripts = deque(sorted(transcripts, reverse=True,
								key=lambda t: (transcripts[t].ba_total_coverage,
											   transcripts[t].ba_coverage_fold)))
	chosen = 0
	while (sorted_transcripts):
		current = sorted_transcripts.popleft()
		for aln in transcripts[current].get_alignments():
			if aln.is_best and not aln.read_assigned:
				aln.read_assigned = True
				aln.is_chosen = True
				chosen += 1
				for per_read_aln in read_aln[aln.read_id]:
					per_read_aln.read_assigned = True
		transcripts[current].join_assigned_regions()
	return dict(filter(lambda item: item[1].get_assigned_read_count(), transcripts.iteritems()))

def remove_orthologue_strains(species2transcript):
	for spec, transcripts in species2transcript.iteritems():
		strain2transcript = defaultdict(list)
		for t in transcripts:
			strain2transcript[t.tax_id].append(t)
		# if only single strain is present, just continue
		if len(strain2transcript) == 1:
			continue
		strains = sorted(strain2transcript, key=lambda s: len(strain2transcript[s]), reverse=True)
		ref_strain = strains[0]
		for strain in strains[1:]:
			transcript2orthologue = dict()
			for transcript in strain2transcript[strain]:
				orthologue = None
				max_similarity = 0.
				for ref_transcript in strain2transcript[ref_strain]:
					similarity = simple_cds_orthologue_test(transcript, ref_transcript)
					if similarity > 0.5:
						if similarity > max_similarity:
							max_similarity = similarity
							orthologue = ref_transcript
				if orthologue is not None:
					transcript2orthologue[transcript] = orthologue
			if len(transcript2orthologue) == len(strain2transcript[strain]):
				species2transcript[spec] = filter(lambda cds: cds.tax_id != strain, species2transcript[spec])
	return species2transcript

def find_hypothetical(transcripts):
	hypothetical = map(lambda t: t.id(),
					   filter(lambda t: t.product == 'hypothetical_protein',
					   		  transcripts.itervalues()))
	return hypothetical

def filter_by_coverage_fold(transcripts, coverage_threshold, fold_threshold):
	new_transcripts = filter(lambda t: t.assigned_total_coverage >= coverage_threshold,
									filter(lambda t: t.assigned_coverage_fold >= fold_threshold, transcripts.itervalues()))
	return dict(map(lambda t: (t.get_id(), t), new_transcripts))

def lca_read_assign(transcripts, read_aln, tax_tree, max_alns):
	total_reads = len(read_aln)
	tax2reads = defaultdict(list)

	for read, alns in read_aln.iteritems():
		# TODO: This is BLAST specific, should not be!!
		#
		# First filter out those whose bitscore is below 50.
		# new_alns = filter(lambda aln: aln.bitscore > 50., alns)
		# Next, let's filter out those whose e-value is above 0.01
		# new_alns = filter(lambda aln: aln.e_value < 0.01, new_alns)
		# ten_perc_thresh = new_alns[0].bitscore * 0.9
		# new_alns = filter(lambda aln: aln.bitscore > ten_perc_thresh, new_alns)
		aln_count = len(alns)
		alns_to_take = min(aln_count, max_alns)
		new_alns = alns[:alns_to_take]
		taxa = map(lambda aln: transcripts[aln.target_id].tax_id, new_alns)
		lca_tax = tax_tree.find_lca(taxa)
		tax2reads[lca_tax].append(read)

	min_read_count = 0.001 * total_reads
	# Let's try this: we find which nodes do not have enough reads
	# assigned to them, then let's try to assign them to a higher
	# ranking node.
	# Step One: Separate the ones who have and the ones that do not
	# have enough reads
	underrepresented = set(filter(lambda tax: len(tax2reads[tax]) < min_read_count, tax2reads))
	reported = set(tax2reads) - underrepresented

	under_tmp = set(underrepresented)
	for tax in underrepresented:
		# This means we've already analyzed it.
		if tax not in under_tmp:
			continue
		# let's ignore the root and the 'unindentified' nodes
		if tax in (0, 1, -1):
			continue
		# Next, let's check if this guy has any children
		# in the underrepresented pack.
		if tax in tax_tree.child_nodes:
			children = set(tax_tree.get_all_children(tax)) & set(under_tmp)
			for child_tax in children:
				tax2reads[tax].extend(tax2reads[child_tax])
				under_tmp.remove(child_tax)
		# If the node has enough reads now, let's move on
		if len(tax2reads[tax]) > min_read_count:
			reported.add(tax)
			under_tmp.remove(tax)
			continue

		# If not, well, let's find the boy's ancestors.
		# First case: At least one of his ancestors is reported.
		found = False
		lineage = list(tax_tree.get_lineage(tax)).reverse()
		# Has no parents. Just throw him out and move on.
		if not lineage:
			under_tmp.remove(tax)
			continue

		for node in lineage[1:]:
			if node in reported:
				tax2reads[node].extend(tax2reads[tax])
				under_tmp.remove(tax)
				found = True
				break
		if not found:
			continue

	return dict(map(lambda tax: (tax, tax2reads[tax]), reported))
