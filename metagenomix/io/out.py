import os
import json
import shutil
from collections import namedtuple, defaultdict, OrderedDict
from contextlib import nested

from Bio import SeqIO
import metagenomix.utils as utils
from metagenomix.io.db import get_cds_id

_tax_attr = ('tax', 'organism_name', 'read_count', 'read_perc')
class _TaxInfo(namedtuple('TaxInfo', _tax_attr)):

	def __iter__(self):
		for attr in _tax_attr:
			yield attr, getattr(self, attr)

	def __dict__(self):
		return {k:v for k, v in iter(self)}

	def __str__(self):
		return ','.join([str(self.tax), self.organism_name,
						 str(self.read_count), str(self.read_perc)])


'''Iterate over tax assignment info'''
def _iter_tax_info(tax2reads, tax_tree, output_reads):
	total_reads = float(sum(map(lambda l: len(l), tax2reads.itervalues())))
	for tax, reads in tax2reads.iteritems():
		count = len(reads)
		perc = count / total_reads
		org_name = tax_tree.get_org_name(tax)
		ti = _TaxInfo(tax, org_name, count, perc)
		yield ti


'''Create and store json transcript information file for each identified tax
Input:
	species2transcript: maps species tax to found transcripts
	output_directory: where to store the files
Output:
	List of created files
'''
def json_transcripts(species2transcript, output_directory, tax_tree):
	files = []
	for species, transcripts in species2transcript.iteritems():
		covfold = OrderedDict({'hypothetical_x': [], 'hypothetical': [],
							   'non_hypothetical_x': [], 'non_hypothetical': []})
		names = defaultdict(list)
		org_name = utils.get_valid_filename(tax_tree.get_org_name(species))
		for t in transcripts:
			if t.product == 'hypothetical_protein':
				covfold['hypothetical_x'].append(t.assigned_total_coverage)
				covfold['hypothetical'].append(t.assigned_coverage_fold)
			else:
				covfold['non_hypothetical_x'].append(t.assigned_total_coverage)
				covfold['non_hypothetical'].append(t.assigned_coverage_fold)
				names['gene'].append(t.gene)
				names['product'].append(t.product)

		covfold_json = '{0}{1}{2}-covfold.json'.format(output_directory, os.path.sep, org_name)
		names_json = '{0}{1}{2}-names.json'.format(output_directory, os.path.sep, org_name)
		files.append((covfold_json, names_json))
		with open(covfold_json, 'w') as fout:
			fout.write(json.dumps(covfold, indent=4))
		with open(names_json, 'w') as fout:
			fout.write(json.dumps(names, indent=4))
	return files

'''Output gene expression data for all species

The data is stored under gene_expression directory of the output.
For each tax of the identified species a separate directory is created
and is filled with strain-specific gene expression data.
'''
def output_species_transcripts(species2transcript, output_directory, tax_tree, assigned=False):
	if not species2transcript:
		return
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)
	if assigned:
		get_coverage = lambda t: t.assigned_total_coverage
		get_fold = lambda t: t.assigned_coverage_fold
		get_read_count = lambda t: t.get_assigned_read_count()
	else:
		get_coverage = lambda t: t.total_coverage
		get_fold = lambda t: t.coverage_fold
		get_read_count = lambda t: t.get_num_of_aln_reads()

	key_attribs = species2transcript.itervalues().next()[0].get_key_attributes(assigned)
	cvs_header = ','.join(attr.upper() for attr in key_attribs)

	i = 0
	for spec_tax, transcripts in sorted(species2transcript.items(), reverse=True, key=lambda i: len(i[1])):
		if spec_tax in (0, -1):
			print 'Warning, taxID in (0, -1), this must be corrected!'
			continue
		i += 1
		strain2transcript = defaultdict(list)
		count = len(transcripts)

		fname = utils.get_valid_filename(tax_tree.get_org_name(spec_tax))
		species_outdir = output_directory + os.path.sep + '%d-%s' % (i, '%s (%d)' % (fname, count))
		if not os.path.exists(species_outdir):
			os.makedirs(species_outdir)

		for t in transcripts:
			strain2transcript[t.tax_id].append(t)

		for strain, strain_trans in strain2transcript.iteritems():
			org_name = tax_tree.get_org_name(strain)
			fname = utils.get_valid_filename(org_name)
			file_name = '%s%s%s.txt' % (species_outdir, os.path.sep, '%s (%d)' % (fname, len(strain_trans)))
			org_name = tax_tree.nodes[strain].organism_name
			sorted_transcripts = sorted(strain_trans, reverse=True,
									key=lambda t: (get_coverage(t), get_fold(t)))
			with open(file_name, 'w') as fout:
				fout.write('#species %s\n\n' % org_name)
				fout.write('#%s\n' % cvs_header)
				for trans in sorted_transcripts:
					fout.write('%s\n' % trans.get_line_descr(key_attribs, ','))

''' General species-gene expression output of metagenomix

The information retrieved via metagenomix is stored via this function
in json format. It marks all taxa identified, their corresponding
transcripts and all the reads mapped to each of the transcripts.
'''
def spec_transcript_read_json(species2transcript, tax_tree, json_output):
	spec_list = []
	for spec, transcripts in species2transcript.iteritems():
		transcript_list = []
		spec_list.append({'taxon': spec,
						  'organism': tax_tree.get_org_name(spec),
						  'rank': tax_tree.get_org_rank(spec),
						  'transcripts': transcript_list})
		for transcript in transcripts:
			transcript_list.append({'id': transcript.get_id(),
									'product': transcript.product,
									'reads': transcript.get_assigned_reads()})
	with open(json_output, 'w') as fout:
		fout.write(json.dumps(spec_list, indent=4))

''' Dumps tax info for all taxa identified in json format '''
def write_t2r_json(tax2reads, tax_tree, fname, output_reads):
	taxa = []
	for taxinfo in _iter_tax_info(tax2reads, tax_tree, output_reads):
		taxa.append(dict(taxinfo))
	with open(fname, 'w') as fout:
		fout.write(json.dumps(taxa, indent=4, sort_keys=True))

''' Dumps tax info for all taxa identified in csv format '''
def write_t2r_csv(tax2reads, tax_tree, fname, output_reads):
	with open(fname, 'w') as fout:
		for taxinfo in _iter_tax_info(tax2reads, tax_tree):
			fout.write('%s\n' % str(taxinfo))

''' Output tax info for specified format '''
def output_tax2reads(tax2reads, tax_tree, fname, out_format, output_reads=False):
	if out_format not in ('json', 'csv'):
		raise ValueError('%s not a valid tax2reads format')
	if out_format == 'json':
		write_t2r_json(tax2reads, tax_tree, fname, output_reads)
	elif out_format == 'csv':
		write_t2r_csv(tax2reads, tax_tree, fname, output_reads)

''' Outputs assigned read count for all identified taxa '''
def species_readcount(species2transcript, tax_tree, fname):
	species2readcount = {}
	for spec, transcripts in species2transcript.iteritems():
		if spec in (0, -1):
			org_name = 'unassigned'
		else:
			org_name = tax_tree.nodes[spec].organism_name
		read_count = reduce(lambda x, y: x + y, map(lambda t: t.get_assigned_read_count(), transcripts))
		species2readcount[org_name] = read_count
	with open(fname, 'w') as fout:
		fout.write('\n'.join(['%s,%d' % (spec, count) for spec, count in species2readcount.iteritems()]))

''' Outputs a bag of read identifiers to newline separated file '''
def output_reads(reads, out_dir, fname, sep='\n'):
	out_file = out_dir + os.path.sep + fname
	with open(out_file, 'w') as fout:
		fout.write(sep.join(reads))

def extract_reads(reads, original_fasta, output_file, file_type):
	reads = set(reads)
	with nested(open(original_fasta), open(output_file, 'w')) as (fin, fout):
		records = SeqIO.parse(fin, file_type)
		for rec in records:
			if rec.name in reads:
				SeqIO.write(rec, fout, file_type)

''' Outputs the reconstructed portion of a taxonomic tree in json format

The structure follows the principles of the taxonomic tree. The root node
is the ancestor of all other nodes, and taxonomic parent is ancestor of
its taxonomic children.
'''

def tax_tree_from_clusters(clusters, tax_tree, fname):
	tax2info = dict()
	all_nodes = set()
	tax2info['organism_name'] = 'root'
	tax2info['tax'] = 1
	tax2info['children'] = []
	tax2info['count'] = 0.
	taxa = set([c.lca for c in clusters])
	tax2cluster = {c.lca: c for c in clusters}
	total_reads = float(sum([len(c.ureads) for c in clusters]))
	for cluster in clusters:
		tax = cluster.lca
		if tax in (0, -1):
			continue
		all_nodes.update(tax_tree.get_lineage(tax))

	def _update_node(root):
		root_tax = root['tax']
		children = filter(lambda node: node in tax_tree.child_nodes[root_tax], all_nodes)
		for child in children:
			if child in taxa:
				count = len(tax2cluster[child].ureads) / total_reads * 100.
				count = round(count, 2)
			else:
				count = 0.
			new_dict = {'tax': child,
						'organism_name': tax_tree.nodes[child].organism_name,
						'children': [],
						'count': count}
			root['children'].append(_update_node(new_dict))
		return root

	def _update_counts(root):
		if not tax_tree.child_nodes[root['tax']]:
			return root['count']
		else:
			root['count'] += round(sum(map(lambda c: _update_counts(c), root['children'])), 2)
			return root['count']

	root = _update_node(tax2info)
	_update_counts(root)

	with open(fname, 'w') as fout:
		fout.write(json.dumps(tax2info, indent=4, sort_keys=True))

def output_tax_tree(species2transcript, tax_tree, fname):
	tax2info = dict()
	all_nodes = set()
	tax2info['organism_name'] = 'root'
	tax2info['tax'] = 1
	tax2info['children'] = []
	tax2info['count'] = 0
	for spec_tax in species2transcript.iterkeys():
		if spec_tax in (0, -1):
			continue
		all_nodes.update(tax_tree.get_lineage(spec_tax))
	def _update_node(root):
		root_tax = root['tax']
		children = filter(lambda node: node in tax_tree.child_nodes[root_tax], all_nodes)
		for child in children:
			if child in species2transcript:
				count = sum(map(lambda t: len(t.get_alignments()), species2transcript[child]))
			else:
				count = 0
			new_dict = {'tax': child,
						'organism_name': tax_tree.nodes[child].organism_name,
						'children': [],
						'count': count}
			root['children'].append(_update_node(new_dict))
		return root

	def _update_counts(root):
		if root['count'] != 0:
			return root['count']
		else:
			root['count'] += sum(map(lambda c: _update_counts(c), root['children']))
			return root['count']

	root = _update_node(tax2info)
	_update_counts(root)

	with open(fname, 'w') as fout:
		fout.write(json.dumps(tax2info, indent=4, sort_keys=True))

''' Fills a d3 template with data and writes to disk '''
def default_d3_plot(input_json, template_html, output_html, filename_repetitions=1):
	fname  = input_json.split(os.path.sep)[-1]
	args = (fname,) * filename_repetitions
	with open(template_html, 'r') as fin:
		html = fin.read()
		html = html % args
	with open(output_html, 'w') as fout:
		fout.write(html)
	out_dir = os.path.sep.join(output_html.split(os.path.sep)[:-1])
	src = input_json
	dst = '%s%s%s' % (out_dir, os.path.sep, fname)
	if src != dst:
		shutil.copyfile(src, dst)

''' Creates a c3js scatter plot from supplied data '''
def c3_scatter(json_covfilt, json_names, template_html, output_html):
	with open(template_html) as fin:
		html = fin.read()
	with open(json_covfilt) as fin:
		covfilt = fin.read()
	with open(json_names) as fin:
		names = fin.read()

	html = html % (names, covfilt)
	with open(output_html, 'w') as fout:
		fout.write(html)

def _get_reads(sp_data):
	transcripts = sp_data['transcripts']
	reads = reduce(lambda a, b: a + b, map(lambda t: t['reads'], transcripts))
	return reads

def _get_products(sp_data):
	transcripts = sp_data['transcripts']
	products = map(lambda t: t['product'], transcripts)
	return products

def _get_species_by_tax(species, tax):
	for spec in species:
		if spec['taxon'] == tax:
			return spec

''' Writes to json file a result of metagenomix-metasim read assignment comparison '''
def read_comparison_json(metagenomix_json, metasim_json, output_json, tax_tree):
	with open(metagenomix_json) as fin:
		metagenomix_data = json.loads(fin.read())
	with open(metasim_json) as fin:
		metasim_data = json.loads(fin.read())
	species_list = []
	metagenomix_species = set(map(lambda d: d['taxon'], metagenomix_data))
	metasim_species = set(map(lambda d: d['taxon'], metasim_data))
	all_species = metasim_species | metagenomix_species

	for species in all_species:
		org_name = tax_tree.nodes[species].organism_name
		if species in metasim_species:
			msim_data = _get_species_by_tax(metasim_data, species)
			msim_reads = set(_get_reads(msim_data))
		else:
			msim_reads = set()

		if species in metagenomix_species:
			pipe_data = _get_species_by_tax(metagenomix_data, species)
			pipe_reads = set(_get_reads(pipe_data))
		else:
			pipe_reads = set()

		total = len(msim_reads)
		true_pos = len(msim_reads & pipe_reads)
		true_neg = len(msim_reads - pipe_reads)
		false_pos = len(pipe_reads - msim_reads)
		species_list.append({'total': total,
						  'true-positive': true_pos,
						  'false-negative': true_neg,
						  'false-positive': false_pos,
						  'species': org_name})
	with open(output_json, 'w') as fout:
		fout.write(json.dumps(species_list))

''' Writes to json file a result of metagenomix-metasim transcript composition comparison '''
def transcript_comparison_json(metagenomix_json, metasim_json, output_json, tax_tree):
	with open(metagenomix_json) as fin:
		metagenomix_data = json.loads(fin.read())
	with open(metasim_json) as fin:
		metasim_data = json.loads(fin.read())
	product_list = []
	metagenomix_species = set(map(lambda d: d['taxon'], metagenomix_data))
	metasim_species = set(map(lambda d: d['taxon'], metasim_data))
	all_species = metasim_species | metagenomix_species

	for species in all_species:
		org_name = tax_tree.nodes[species].organism_name
		if species in metasim_species:
			msim_data = _get_species_by_tax(metasim_data, species)
			msim_products = set(_get_products(msim_data))
		else:
			msim_products = set()

		if species in metagenomix_species:
			pipe_data = _get_species_by_tax(metagenomix_data, species)
			pipe_products = set(_get_products(pipe_data))
		else:
			pipe_products = set()

		total = len(msim_products)
		true_pos = len(msim_products & pipe_products)
		true_neg = len(msim_products - pipe_products)
		false_pos = len(pipe_products - msim_products)
		product_list.append({'total': total,
							 'true-positive': true_pos,
							 'false-negative': true_neg,
							 'false-positive': false_pos,
							 'species': org_name})
		with open(output_json, 'w') as fout:
			fout.write(json.dumps(product_list))

def output_clusters(clusters, output_file):
	with open(output_file, 'w') as fout:
		for cluster in clusters:
			fout.write('{}\n'.format(str(cluster)))