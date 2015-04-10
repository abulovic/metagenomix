import re
import json
from collections import defaultdict

from Bio import SeqIO

from metagenomix.io.db import parse_cds_header, get_cds_id

def get_read_count(metasim_fasta, tax_tree):
	'''Returns dictionary species_tax: read_count'''
	species2readcount = defaultdict(int)
	with open(metasim_fasta) as fin:
		records = SeqIO.parse(fin, 'fasta')
		for record in records:
			data = record.description.split('|')
			tax_idx = data.index('taxon') + 1
			tax = int(data[tax_idx])
			species_tax = tax_tree.get_parent_with_rank(tax, 'species')
			if species_tax in (0, -1):
				org_name = 'unassigned'
			else:
				org_name = tax_tree.nodes[species_tax].organism_name
			species2readcount[org_name] += 1
	return species2readcount

def synth_dataset_profile(metasim_fasta, tax_tree, json_output):
	'''Creates a JSON file	of format:
	{species_tax:
		{GeneID1: [
			readID1,
			readID2
			]
		 GeneID2: [
		 	readID35,
		 	readID38
		 	]
		 }
	}'''
	spec2transcripts = defaultdict(lambda: defaultdict(list))
	id2product = {}
	with open(metasim_fasta) as fin:
		records = SeqIO.parse(fin, 'fasta')
		for record in records:
			try:
				header = re.findall(r'"(.*?)"', record.description)[0]
			except Exception, e:
				import pdb
				pdb.set_trace()
			data = parse_cds_header(header)
			tax = int(data['taxon'])
			species_tax = tax_tree.get_parent_with_rank(tax, 'species')
			if species_tax in (0, -1):
				continue
			_id = get_cds_id(data)
			id2product[_id] = data['product']
			read_id = record.name
			spec2transcripts[species_tax][_id].append(read_id)
	species_list = []
	for spec_tax, transcripts in spec2transcripts.iteritems():
		transcript_list = []
		_data = {'taxon': spec_tax, 'transcripts': transcript_list}
		for tid, reads in transcripts.iteritems():
			product = id2product[tid]
			transcript_list.append({'id': tid,
									'product': product,
									'reads': reads})
		species_list.append(_data)

	with open(json_output, 'w') as fout:
		fout.write(json.dumps(species_list, indent=4))

