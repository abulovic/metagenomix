'''
A module offering functionality to parse the megablast output.
'''
import re
from collections import defaultdict
import lxml.etree as ET

from Bio.Blast import NCBIXML
from metagenomix.sequence.seq import CDS, BlastAlignment, TargetSeq
from metagenomix.utils import timeit, oneline_timeit, progressbar
from metagenomix.io.db import parse_cds_header
from metagenomix.sequence.seq import create_sequence, get_seq_id

_location_pattern = re.compile("(\d+)")

def _parse_n_lines(infile, n_lines, read_alns, target_seqs, db_type):

	new_reads = set()
	for i, line in enumerate(infile):
		if i >= n_lines:
			return new_reads
		if line.startswith('#'):
				continue
		data = line.strip().split('\t')
		read_id, target, _, _, _, _, qstart, qend, tstart, tend, e_value, bitscore = data
		if read_id not in read_alns:
			new_reads.add(read_id)
		tstart, tend = int(tstart), int(tend)
		if tstart > tend:
			tstart, tend = tend, tstart
		if db_type == 'cds':
			header_data = parse_cds_header(target)
			locations = re.findall(_location_pattern, header_data['location'])
			start, stop = locations[0], locations[-1]
			length = abs(int(start) - int(stop))
		else:
			length = -1

		target_seq = create_sequence(target, length, db_type, detailed=False)
		target_id = target_seq.get_id()
		if target_id not in target_seqs:
			target_seqs[target_id] = target_seq
		target_seq = target_seqs[target_id]

		if target_id not in [aln.target_id for aln in read_alns[read_id]]:
			aln = BlastAlignment(read_id, target_id, qstart, qend, tstart, tend, e_value, bitscore)
			read_alns[read_id].append(aln)
			#target_seq.alignment.add_location(read_id, tstart, tend)

	return new_reads

def _filter_low_scoring(read_alns, bitscore, evalue, new_reads):
	total = 0
	filtered = 0
	for read in new_reads:
		alns = read_alns[read]
		total += len(alns)
		alns[:] = filter(lambda aln: aln.bitscore > bitscore, alns)
		alns[:] = filter(lambda aln: aln.e_value < evalue, alns)
		filtered += len(alns)

def annotate_targets(read_alns, target_seqs):
	with timeit('Removing reads with 0 alignments'):
		zero_aln_reads = filter(lambda r: len(read_alns[r])==0, read_alns)
		for read in zero_aln_reads:
			read_alns.pop(read)
		print 'Removed %d reads' % len(zero_aln_reads)

	# Annotate the best alignment for all the reads.
	with timeit('Adding alignments to target sequences'):
		total_reads = len(read_alns)
		step = max(1, total_reads/100)
		print 'Number of reads to annotate:', total_reads
		progressbar(0, start=True)
		for i, (read, alns) in enumerate(read_alns.iteritems()):
			if i % step == 0:
				progressbar(i/step)
			alns.sort(reverse=True)
			max_score = alns[0].bitscore
			for aln in alns[:max(1, int(0.1*len(alns)))]:
				if aln.bitscore == max_score:
					aln.is_best = True
				else:
					aln.is_best = False
		progressbar(100, end=True)

	with timeit('Functional transcript annotation'):
		for read, alns in read_alns.iteritems():
			for aln in alns:
				target = target_seqs[aln.target_id]
				if aln.is_best:
					target.num_of_ba_reads += 1
				target.num_of_aln_reads += 1
				target.add_alignment(aln)

	with timeit('Removing 0 alignments target seqs'):
		zero_aln_targets = filter(lambda t: len(target_seqs[t].get_alignments()) == 0, target_seqs.iterkeys())
		for target in zero_aln_targets:
			target_seqs.pop(target)
		print 'Removed %d targets.' % len(zero_aln_targets)

	for target in target_seqs.itervalues():
		target.join_regions()
		target.join_ba_regions()

	return read_alns, target_seqs


def parse_tab_delimited(input_file, db_type, at_once=1e5, detailed=False,
						read_alns=defaultdict(list), target_seqs={},
						filter_low_scoring=True, annotate=True, entry_cnt=None):

	if entry_cnt is None:
		with timeit('Line count'):
			with open(input_file) as f:
				line_count = sum(1 for line in f)
	else:
		line_count = entry_cnt
	rolls = int(line_count / at_once + 1)

	with timeit('Pure file parsing'):
		with open(input_file, 'r') as fin:
			for i in xrange(rolls):
				with oneline_timeit('Step %d/%d' % (i, rolls)):
					new_reads = _parse_n_lines(fin, at_once, read_alns, target_seqs, db_type)
					if filter_low_scoring:
						_filter_low_scoring(read_alns, 100., 0.01, new_reads)

	if not annotate:
		return read_alns, target_seqs
	else:
		return annotate_targets(read_alns, target_seqs)

def parse_xml(input_file, db_type, annotate=True, detailed=False, entry_cnt=None):
	target_seqs = {}
	read_alns = defaultdict(list)
	if entry_cnt is None:
		entry_cnt = 0
		with open(input_file) as fin:
			for line in fin:
				if line.strip() == "<Iteration>":
					entry_cnt += 1

	step = max(1, entry_cnt / 100)

	with open(input_file) as fin:
		context = ET.iterparse(fin, events=('end',), tag='Iteration')
		i = 0
		progressbar(0, start=True)
		for event, elem in context:
			i += 1
			if i % step == 0:
				progressbar(i/step)
			iter_children = {e.tag: e for e in elem.getchildren()}
			# Blast outputs a message if it has no alns
			if 'Iteration_message' in iter_children:
				continue
			else:
				hits = iter_children['Iteration_hits'].getchildren()
				read_id = iter_children['Iteration_query-def'].text
				for hit in hits:
					hit_children = {e.tag: e for e in hit.getchildren()}
					target_header = hit_children['Hit_id'].text
					#target_header = hit_children['Hit_def'].text
					tid = get_seq_id(target_header, db_type)
					if tid not in target_seqs:
						target_len = int(hit_children['Hit_len'].text)
						target_seq = create_sequence(target_header, target_len, db_type, detailed)
						target_seqs[tid] = target_seq
					else:
						target_seq = target_seqs[tid]

					hsp = hit_children['Hit_hsps'].getchildren()[0]
					hsp_data = {e.tag: e for e in hsp.getchildren()}
					tstart = int(hsp_data['Hsp_hit-from'].text)
					tend = int(hsp_data['Hsp_hit-to'].text)
					qstart = int(hsp_data['Hsp_query-from'].text)
					qend = int(hsp_data['Hsp_query-to'].text)
					evalue = float(hsp_data['Hsp_evalue'].text)
					bitscore = float(hsp_data['Hsp_bit-score'].text)
					query = hsp_data['Hsp_qseq'].text
					target = hsp_data['Hsp_hseq'].text

					aln = BlastAlignment(read_id, tid, qstart, qend, tstart, tend, evalue, bitscore)
					read_alns[read_id].append(aln)

					if detailed:
						target_seq.alignment.add_location(read_id, tstart, tend, query, target)
					else:
						target_seq.alignment.add_location(read_id, tstart, tend)
			elem.clear()
			for ancestor in elem.xpath('ancestor-or-self::*'):
				while ancestor.getprevious() is not None:
					del ancestor.getparent()[0]
		progressbar(100, end=True)
		del context

	return annotate_targets(read_alns, target_seqs)

def get_entry_cnt_xml(aln_file):
	reads = 0
	with open(aln_file) as fin:
		for l in fin:
			if l.strip() == '<Iteration>':
				reads += 1
	return reads

def get_entry_cnt_tab(aln_file):
	reads = 0
	with open(aln_file) as fin:
		for l in fin:
			reads += 1
	return reads

def extract_unique_reads(input_file):
	reads = set()
	with open(input_file) as fin:
		for line in fin:
			if line.startswith('#'):
				continue
			read_id = line.split()[0]
			reads.add(read_id)
	return reads
