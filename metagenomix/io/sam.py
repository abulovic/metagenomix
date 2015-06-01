import os
import time
from collections import defaultdict

import pysam
import numpy as np

from metagenomix.sequence.seq import CDS, SamAlignment, TargetSeq, Genome
from metagenomix.utils import timeit, get_appropriate_file_size
from metagenomix.io.db import parse_cds_header, parse_genome_header
from metagenomix.sequence.seq import create_sequence
from metagenomix.biodb.data_access import DataAccess

def annotate_targets(read_alns, target_seqs):

	with timeit('Semantic target sequence annotation'):
		for read, alns in read_alns.iteritems():
			for aln in alns:
				target = target_seqs[aln.target_id]
				if aln.is_best:
					target.num_of_ba_reads += 1
				target.num_of_aln_reads += 1
				target.add_alignment(aln)

		for target in target_seqs.itervalues():
			target.join_regions()
			target.join_ba_regions()

	return read_alns, target_seqs


def parse_cds_sam(input_file, db_type, binary=False, annotate=True,
				  read_alns=defaultdict(list), target_seqs={}, entry_cnt=None, detailed=False):
	with timeit('Pysam load of the %s file (%s)' % (input_file.split(os.path.sep)[-1],
													get_appropriate_file_size(input_file))):
		if binary:
			samfile = pysam.Samfile(input_file, 'rb')
		else:
			samfile = pysam.Samfile(input_file, 'r')

	if entry_cnt is not None:
		_1perc_reads = max(1, entry_cnt / 100)
	data_access = DataAccess()
	unmapped = set()
	lengths = np.array(samfile.lengths)
	with timeit('Sequence iteration'):
		for i, ar in enumerate(samfile.fetch()):
			if entry_cnt is not None:
				if i % _1perc_reads == 0:
					print '%d %%' % (i / _1perc_reads)
			else:
				if i % 10000 == 0:
					print 'Processed %d aligned reads.' % i
			if ar.is_unmapped:
				unmapped.add(ar.qname)
				continue

			read_id = ar.qname.split()[0]
			tstart, tend = ar.aend - ar.alen, ar.aend
			try:
				qstart, qend = ar.qstart, ar.qend
			except SystemError as e:
				# pysam throws a System Error! :(
				continue
			if ar.tid == -1:
				continue
			target_name = samfile.getrname(ar.tid)
			target_len = lengths[ar.tid]
			target_seq = create_sequence(target_name, target_len, db_type)
			if target_seq.get_id() not in target_seqs:
				target_seqs[target_seq.get_id()] = target_seq
			else:
				target_seq = target_seqs[target_seq.get_id()]

			target_seq.alignment.add_location(read_id, tstart, tend)
			aln = SamAlignment(read_id, target_seq.get_id(), qstart, qend, tstart,
							   tend, ar.mapq)
			if not ar.is_secondary:
				aln.is_best = True
			read_alns[read_id].append(aln)

	if not annotate:
		return read_alns, target_seqs
	else:
		return annotate_targets(read_alns, target_seqs)

	samfile.close()


def get_entry_cnt_sam(aln_file):
	entries = 0
	with open(aln_file) as fin:
		for l in fin:
			if not l.startswith("@SQ"):
				entries += 1
	return entries

def extract_unique_reads(input_file):
	if input_file.endswith('sam'):
		samfile = pysam.Samfile(input_file, 'r')
	elif input_file.endswith('bam'):
		samfile = pysam.Samfile(input_file, 'rb')

	reads = set()
	for record in samfile.fetch():
		if record.is_unmapped:
			continue
		read_id = record.qname.split()[0]
		reads.add(read_id)
	return reads
