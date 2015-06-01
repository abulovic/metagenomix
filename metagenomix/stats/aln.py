import os
from collections import Counter, defaultdict

import numpy as np

from metagenomix.stats.otu import get_rank_distribution


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

