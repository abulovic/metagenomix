import os
import json
from collections import defaultdict

_invalid_taxa = (0, -1)


def get_rank_distribution(species2transcript, tax_tree, rank):
	rank2count = defaultdict(int)
	if isinstance(species2transcript, list):
		for cluster in species2transcript:
			if tax_tree.get_org_rank(cluster.lca) == rank:
				parent = cluster.lca
			else:
				parent = tax_tree.get_parent_with_rank(cluster.lca, rank)
			rank2count[parent] += len(cluster.ureads)
	else:
		for tax, transcripts in species2transcript.iteritems():
			parent = tax_tree.get_parent_with_rank(tax, rank)
			rank2count[parent] += len(transcripts)
	return rank2count

def _write_tax_distr_line(fout, tax, count, total, tax_tree):
	if tax in _invalid_taxa:
		name = 'unassigned'
	else:
		name = tax_tree.nodes[tax].organism_name
	fout.write('%d,%s,%d,%.3f\n' % (tax, name, count,
									(float(count)/total*100.)))

def simple_taxa_transcript_distribution(species2transcript, output_dir,
										tax_tree, naming_prefix, additional_ranks=('phylum',)):

	if naming_prefix == 'clusters':
		spec2count = {c.lca: len(c.ureads) for c in species2transcript}
	else:
		spec2count = dict(map(lambda tax: (tax, len(species2transcript[tax])),
							  species2transcript.keys()))
	newrank2count = defaultdict(int)
	total_transcripts = sum(spec2count.values())
	csv_files = []
	json_files = []

	file_name = output_dir + os.path.sep + '%s-simple-species-distribution.csv' % naming_prefix
	csv_files.append(file_name)
	with open(file_name, 'w') as fout:
		fout.write('TAX,SPECIES,TRANSCRIPT COUNT,PERCENTAGE IN SAMPLE\n')
		for tax, count in sorted(spec2count.items(), key=lambda i: i[1], reverse=True):
			_write_tax_distr_line(fout, tax, count, total_transcripts, tax_tree)

	for rank in additional_ranks:
		fname = '{0}{1}{2}-{3}-distribution.csv'.format(output_dir, os.path.sep,
														naming_prefix, rank)
		csv_files.append(fname)
		rank2count = get_rank_distribution(species2transcript, tax_tree, rank)
		with open(fname, 'w') as fout:
			fout.write('TAX,SPECIES,TRANSCRIPT COUNT,PERCENTAGE IN SAMPLE\n')
			for tax, count in sorted(rank2count.items(), key=lambda i: i[1], reverse=True):
				_write_tax_distr_line(fout, tax, count, total_transcripts, tax_tree)
		fname = '{0}{1}{2}-{3}-distribution.json'.format(output_dir, os.path.sep,
														naming_prefix, rank)
		json_files.append(fname)
		name2count = {}
		for rank, count in rank2count.iteritems():
			if rank in (0, -1):
				org_name = 'unassigned'
			else:
				org_name = tax_tree.nodes[rank].organism_name
			name2count[org_name] = count
		with open(fname, 'w') as fout:
			fout.write(json.dumps(name2count))
	return csv_files, json_files
