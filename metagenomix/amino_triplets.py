import os
from subprocess import call, Popen, PIPE
from itertools import product
from collections import defaultdict
from random import randint, shuffle

import numpy as np
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACProtein as ProteinAlphabet
from gensim import corpora, models, similarities

from metagenomix.tax import TaxTree, ranks
from metagenomix.utils import time_exec, timeit
from metagenomix.biodb.data_access import DataAccess

def iter_over_species_trips(seq_dir):
	with open(os.path.sep.join([seq_dir, 'name2tax.txt'])) as fin:
			name2tax = dict(map(lambda d: (d[1], int(d[0])), [l.strip().split('|') for l in fin]))
	for f in sorted(os.listdir(seq_dir)):
		if f.endswith('txt') and f != 'name2tax.txt':
			org_name = f[:-4]
			full_f = os.path.sep.join([seq_dir, f])
			trips = []
			if org_name not in name2tax:
				print 'MISSING ', org_name
				continue
			with open(full_f) as fin:
				for seq in fin:
					seq = seq.strip()
					trips.extend([seq[x:x+3] for x in xrange(0, len(seq)-3)])
			yield org_name, trips

def iter_over_species_files(seq_dir):
	with open(os.path.sep.join([seq_dir, 'name2tax.txt'])) as fin:
			name2tax = dict(map(lambda d: (d[1], int(d[0])), [l.strip().split('|') for l in fin]))
	for f in sorted(os.listdir(seq_dir)):
		if f.endswith('txt') and f != 'name2tax.txt':
			org_name = f[:-4]
			full_f = os.path.sep.join([seq_dir, f])
			if org_name not in name2tax:
				print 'MISSING ', org_name
				continue
			yield full_f


class SpeciesCorpus(object):
	def __init__(self, seq_dir, freq_tool):
		super(SpeciesCorpus, self).__init__()
		self.seq_dir = seq_dir
		self.freq_tool = freq_tool
		with open(os.path.sep.join([seq_dir, 'name2tax.txt'])) as fin:
			self.name2tax = dict(map(lambda d: (d[1], int(d[0])), [l.strip().split('|') for l in fin]))
		self.index2tax = {}
		with timeit('Initializing dictionary'):
			self.dictionary = corpora.Dictionary((data for data in iter_over_species_files(self.seq_dir)), freq_tool=freq_tool)

	def __iter__(self):
		print 'ITER ITER'
		for fname in iter_over_species_files(self.seq_dir):
		#for i, (org_name, trips) in enumerate(iter_over_species_trips(self.seq_dir)):
		#	self.index2tax[i] = self.name2tax[org_name]
			yield self.dictionary.doc2bow(fname)

	def __next__(self):
		print 'NEXT NEXT'
		for i, (org_name, trips) in enumerate(iter_over_species_trips(self.seq_dir)):
			self.index2tax[i] = self.name2tax[org_name]
			yield self.dictionary.doc2bow(trips)
		raise StopIteration


def close_third_files(files):
	shuffle(files)
	for f in iter(files[0:len(files)/3]):
		if not f.closed:
			f.close()


@time_exec('SVE IKAD')
def nesto():
	out_dir = '/home/abulovic/tmp/'
	org2f = {}
	name2tax = {}
	da = DataAccess()
	old_name = None
	fname = '/home/abulovic/BINNER/database/nr/bacteria.nr_protein.faa'
	with open(fname) as fin:
		recs = SeqIO.parse(fin, 'fasta')
		for rec in recs:
			organism_name = rec.description.split('[')[-1][:-1].replace('/', '_')
			if organism_name not in name2tax:
				gi = int(rec.description.split('|')[1])
				tax = da.get_taxids((gi,), tuple)
				if tax:
					name2tax[organism_name] = tax[0]
			if old_name != organism_name:
				old_name = organism_name
			if organism_name not in org2f:
				try:
					org2f[organism_name] = open(out_dir + organism_name + '.txt', 'w')
				except IOError, e:
					print e, len(filter(lambda f: not f.closed, org2f.values()))
					close_third_files(org2f.values())
					org2f[organism_name] = open(out_dir + organism_name + '.txt', 'w')

			if org2f[organism_name].closed:
				try:
					print e, len(filter(lambda f: not f.closed, org2f.values()))
					org2f[organism_name] = open(out_dir + organism_name + '.txt', 'w+')
				except IOError, e:
					close_third_files(org2f.values())
					org2f[organism_name] = open(out_dir + organism_name + '.txt', 'w+')
			org2f[organism_name].write(str(rec.seq) + '\n')
	for f in org2f.values():
		f.close()
	with open(out_dir + 'name2tax.txt', 'w') as fout:
		for name, tax in name2tax.iteritems():
			fout.write('%d|%s\n' % (tax, name))
