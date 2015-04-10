import os
import sys
import time
import shutil
from contextlib import contextmanager
from collections import defaultdict

from metagenomix.biodb.data_access import DataAccess


class NotSupportedError(Exception):
	pass

@contextmanager
def timeit(func_name):
	print 'Started <%s>.' % func_name
	start =  time.time()
	yield
	stop = time.time()
	print 'Execution of <%s> lasted: %.3f secs.\n' % (func_name, (stop-start))

@contextmanager
def oneline_timeit(message):
	start = time.time()
	yield
	stop = time.time()
	print '%s <%.3f sec>' % (message, (stop - start))

def get_user_reply(message, options, descriptions=None):
	print message
	if descriptions is not None:
		for key in options:
			value = descriptions[key]
			print '  (%s) %s' % (key, value)
	print
	reply = '---'
	options = map(lambda s: s.upper(), options)
	while reply.upper() not in options:
		reply = raw_input('%s > ' % ' / '.join(options))
	return reply.lower()

def table_exists(db, cursor, table_name):
	cursor.execute("SELECT name FROM sqlite_master")
	table_names = map(lambda t: t[0], cursor.fetchall())
	if unicode(table_name) in table_names:
		return True
	else:
		return False

class time_exec(object):
	def __init__(self, message):
		self.message = message

	def __call__(self, func):
		def wrapped(*args):
			print 'Started <%s>.' % self.message
			start = time.time()
			res = func(*args)
			stop = time.time()
			print 'Execution of <%s> lasted: %.3f secs.\n' % (self.message, (stop-start))
			return res
		return wrapped

def progressbar(step, start=False, end=False):
	if not start:
		sys.stdout.write('\r')
	sys.stdout.write('[')
	sys.stdout.write('='*step)
	sys.stdout.write(' '*(100-step))
	sys.stdout.write('] %2d%%' % step)
	if end:
		sys.stdout.write('\n')
	sys.stdout.flush()

class ValueNotSetError(Exception):
	pass

def get_appropriate_file_size(file_name):
	# lets get size in bytes
	size = os.path.getsize(file_name)
	new_size = size / 1024.
	if new_size < 1:
		return "%d B" % size
	size = new_size / 1024.
	if size < 1:
		return "%.2f kB" % new_size
	new_size = size / 1024.
	if new_size < 1:
		return "%.2f MB" % size
	else:
		return "%.2f GB" % new_size

def get_file_type(seq_file):
	seq_file = seq_file.lower()
	ending = seq_file.split('.')[-1]
	if ending in ('fa', 'fasta'):
		return 'fasta'
	elif ending in ('fq', 'fastq'):
		return 'fastq'
	elif ending in ('sam', 'bam'):
		return ending
	elif ending in ('mgb', 'blast', 'blastout', 'm8'):
		return 'blast'
	elif ending == 'xml':
		return 'xml'
	return None

'''
Input:
	tax2reads: species or strain taxa mapped to reads assigned to it
	tax_tree: instance of TaxTree
Output:
	Specified rank taxa mapped to corresponding child reads
'''
def get_rank_read_distribution(tax2reads, tax_tree, rank='species'):
	rank2reads = defaultdict(list)
	for tax, reads in tax2reads.iteritems():
		rank_tax = tax_tree.get_parent_with_rank(tax, rank)
		rank2reads[rank_tax].extend(reads)
	return rank2reads

'''
Input:
	transcripts: species or strain taxa mapped to transcripts assigned to it
	tax_tree: instance of TaxTree
Output:
	Specified rank taxa mapped to corresponding child transcripts
'''
def get_species_transcript_distribution(transcripts, tax_tree):
	species2transcript = defaultdict(list)
	for trans_id, transcript in transcripts.iteritems():
		species_tax = tax_tree.get_parent_with_rank(transcript.tax_id, 'species')
		species2transcript[species_tax].append(transcript)
	return species2transcript


def get_valid_filename(fname):
	new_fname = fname.replace(os.path.sep, '_')
	new_fname = '_'.join(new_fname.split())
	return new_fname

def simple_orthologue_test(g1, g2, l1, l2, p1, p2):
	# first we check if the gene names are the same
	if g1 == g2 and g1 not in ('', 'None') and g1 is not None:
		return 1.
	if not p1 or not p2:
		return 0.
	# if the gene names do not correspond, we check the
	# length ratio and product description overlap
	len_ratio = float(min(l1, l2)) / max(l1, l2)
	if 'hypothetical' in p1 or 'hypothetical' in p2:
		# if either is tagged as a hypothetical protein,
		# the chances for orthologicity are 50-50
		p_product = 0.5
	if not p1 or not p2:
		# if either of them lack the product description,
		# the chances for orthologicity are 50-50
		p_product = 0.5
	else:
		pr1 = set(p1.split('_'))
		pr2 = set(p2.split('_'))
		p_product = len(pr1 & pr2) / float(max(len(pr1), len(pr2)))
	return len_ratio * p_product

def simple_cds_orthologue_test(cds1, cds2):
	return simple_orthologue_test(cds1.gene, cds2.gene,
								  cds1.length, cds2.length,
								  cds1.product, cds2.product)

def retrieve_tax_ids(target_seqs, db_type):
	data_access = DataAccess()
	if db_type == 'cds':
		gi_type = 'nucl_gi'
	else:
		gi_type = 'gi'
	get_gi = lambda t: getattr(t, gi_type)

	with timeit('GI2TAX database querying'):
		gis = set(map(lambda t: get_gi(t), target_seqs.itervalues()))
		if db_type == 'nr':
			tax_ids = data_access.get_taxids(gis, 'prot')
		else:
			tax_ids = data_access.get_taxids(gis, 'nucl')
		missing_gis = set()
		for target in target_seqs.itervalues():
			target.tax_id = tax_ids.get(get_gi(target), -1)
			if target.tax_id == -1:
				missing_gis.add(get_gi(target))

def create_dir(dir_path, erase):
	if os.path.isdir(dir_path):
		if erase:
			shutil.rmtree(dir_path)
			os.makedirs(dir_path)
	else:
		os.makedirs(dir_path)

def get_seq_count(fasta_file):
	line_cnt = 0
	with open(fasta_file) as fin:
		for line in fin:
			if line.startswith('>'):
				line_cnt += 1
	return line_cnt

def reads_from_fasta(fasta_file):
	reads = set()
	with open(fasta_file) as fin:
		for line in fin:
			line = line.strip()
			if not line.startswith('>'):
				continue
			else:
				reads.add(line.split()[0][1:])
	return reads
