import re
from collections import deque

import metagenomix.io.db as db
from metagenomix.utils import ValueNotSetError
from metagenomix.sequence.analysis import SimpleAlignment, DetailedAlignment

_location_pattern = re.compile("(\d+)")


class TargetSeq(object):
	__slots__ = ['accession', 'gi', 'length', 'num_of_aln_reads', 'num_of_ba_reads', 'tax_id',
				 '_preprocessed_regions', '_preprocessed_ba_regions',
				 '_regions_coverage', '_regions_ba_coverage',
				 '_preprocessed_assigned_regions', '_assigned_regions',
				 '_assigned_regions_coverage', 'assigned_read_count',
				 '_alignments', 'regions', 'ba_regions',
				 'covered_regions_length', 'ba_covered_regions_length',
				 'total_coverage', 'coverage_fold',
				 'ba_total_coverage', 'ba_coverage_fold',
				 'assigned_total_coverage', 'assigned_coverage_fold',
				 'assigned_cov_regions_length', 'description',
				 'alignment']

	@classmethod
	def from_header(cls, header, length, detailed=False):
		data = db.parse_nt_header(header)
		return cls(data['accession'], data['gi'], length, data['description'], detailed)

	def get_id(self):
		return self.accession

	def __init__(self, accession, gi, length, description='', tax_id=-1, detailed=False):
		self.accession = accession
		self.gi = int(gi)
		self.length = int(length)
		self.description = description

		self.num_of_aln_reads = 0
		self.num_of_ba_reads = 0
		self._preprocessed_regions = list()
		self._regions_coverage = 0
		self._preprocessed_ba_regions = list()
		self._regions_ba_coverage = 0
		self._preprocessed_assigned_regions = list()
		self._assigned_regions = list()
		self._assigned_regions_coverage = 0
		self._alignments = dict()
		self.assigned_read_count = 0
		self.tax_id = tax_id
		if detailed:
			self.alignment = DetailedAlignment(length)
		else:
			self.alignment = SimpleAlignment(length)

	def get_line_descr(self, keys, separator):
		simple_keys = filter(lambda k: not k.startswith('get'), keys)
		values = [str(getattr(self, key)) for key in simple_keys]
		for k in (set(keys) - set(simple_keys)):
			values.append(str(getattr(self, k)()))
		return separator.join(values)

	def get_key_attributes(self, *args, **kwargs):
		return ('accession', 'gi', 'tax_id', 'length',
				'assigned_total_coverage',
				'assigned_coverage_fold',
				'get_alignments_count',
				'get_best_alignments_count')

	def __str__(self):
		repr_str = []
		repr_str.append('%s (%d)' % (self.get_id(), self.tax_id))
		repr_str.append('Coverage: %.3f' % self.total_coverage)
		repr_str.append('Fold    : %.3f' % self.coverage_fold)
		repr_str.append('Read alignments to CDS: %d' % len(self._preprocessed_regions))
		repr_str.append('Best alignments to CDS: %d' % len(self._preprocessed_ba_regions))
		repr_str.append('Assigned alignments   : %d' % len(self._assigned_regions))
		if len(self._assigned_regions):
			repr_str.append('Assinged coverage     : %.3f' % self.assigned_total_coverage)
			repr_str.append('Assinged coverage fold: %.3f' % self.assigned_coverage_fold)
		return '\n'.join(repr_str)


	def assign_species(self, tax_id):
		self.tax_id = tax_id

	def add_alignment(self, alignment):
		if alignment.read_id not in self._alignments:
			self._alignments[alignment.read_id] = alignment
		self._add_region(alignment.tstart, alignment.tend)
		if alignment.is_best:
			self._add_ba_region(alignment.tstart, alignment.tend)

	def _add_region(self, start, stop):
		self._preprocessed_regions.append((start, stop))
		self._regions_coverage += abs(stop - start)

	def _add_ba_region(self, start, stop):
		self._preprocessed_ba_regions.append((start, stop))
		self._regions_ba_coverage += abs(stop - start)

	def _add_assigned_region(self, start, stop):
		self._preprocessed_assigned_regions.append((start, stop))
		self._assigned_regions_coverage += abs(stop - start)

	def _calculate_stats(self):
		self.covered_regions_length = sum(map(lambda o: o[1] - o[0], self.regions))
		self.coverage_fold = float(self._regions_coverage) / self.covered_regions_length
		self.total_coverage = self.covered_regions_length / float(self.length)

	def _calculate_ba_stats(self):
		self.ba_covered_regions_length = sum(map(lambda o: o[1] - o[0], self.ba_regions))
		if not self.ba_covered_regions_length:
			self.ba_coverage_fold = 0
			self.ba_total_coverage = 0
		else:
			self.ba_coverage_fold = float(self._regions_ba_coverage) / self.ba_covered_regions_length
			self.ba_total_coverage = self.ba_covered_regions_length / float(self.length)

	def _calculate_assigned_stats(self):
		self.assigned_cov_regions_length = sum(map(lambda o: o[1] - o[0], self._assigned_regions))
		if not self.assigned_cov_regions_length:
			self.assigned_coverage_fold = 0
			self.assigned_total_coverage = 0
		else:
			self.assigned_coverage_fold = float(self._assigned_regions_coverage) / self.assigned_cov_regions_length
			self.assigned_total_coverage = self.assigned_cov_regions_length / float(self.length)

	def _get_joint_regions(self, input_regions):
		if len(input_regions) <= 1:
			return input_regions
		pre = list(input_regions)
		pre.sort()
		pre = deque(pre)
		regions = list()
		regions.append(pre.popleft())
		while len(pre):
			r1 = regions.pop()
			r2 = pre.popleft()
			if r1[1] > r2[0]:
				regions.append((r1[0], r2[1]))
			else:
				regions.append(r1)
				regions.append(r2)
		return regions

	def join_regions(self):
		self.regions = self._get_joint_regions(self._preprocessed_regions)
		self._calculate_stats()

	def join_ba_regions(self):
		self.ba_regions = self._get_joint_regions(self._preprocessed_ba_regions)
		self._calculate_ba_stats()

	def join_assigned_regions(self):
		for aln in filter(lambda aln: aln.is_chosen, self.get_alignments()):
			self._add_assigned_region(aln.tstart, aln.tend)
		self._assigned_regions = self._get_joint_regions(self._preprocessed_assigned_regions)
		self._calculate_assigned_stats()

	def get_assigned_read_count(self):
		return len(self._preprocessed_assigned_regions)


	def get_covered_regions_length(self):
		self._return_value('covered_regions_length', -1)

	def _return_value(self, name, alternative):
		if name not in dir(self):
			return alternative
		return getattr(self, name)

	def get_total_coverage(self):
		'''Percentage of nucleotides covered with at least one read'''
		self._return_value('total_coverage', -1)

	def get_tax(self):
		self._return_value('tax', -1)

	def get_num_of_aln_reads(self):
		'''
		Total number of reads aligned to the target, even if it is
		not the best alignment for the read
		'''
		return len(self._preprocessed_regions)

	def get_num_of_ba_reads(self):
		'''
		Number of reads for which this is the best alignment.
		Reads can have more than one 'best' alignment, if their
		alignment score says so.
		'''
		self._return_value('num_of_ba_reads', -1)

	def get_covered_areas(self):
		'''Get the disjunct areas covered by all the aligned reads.'''
		self._return_value('covered_areas', [])

	def get_ba_covered_areas(self):
		'''Get the disjunct areas covered by the best-alignment reads.'''
		self._return_value('ba_covered_areas', [])

	def get_coverage_fold(self):
		'''Get the average nucleotide coverage in the covered areas'''
		self._return_value('coverage_fold', -1)

	def get_ba_coverage_fold(self):
		'''
		Get the average nucleotide coverage in the areas covered
		by the best alignment reads.
		'''
		self._return_value('ba_coverage_fold', -1)

	def get_alignments(self, iterate=False):
		if iterate:
			return self._alignments.itervalues()
		return self._alignments.values()

	def get_assigned_reads(self):
		return self._alignments.keys()

	def get_alignments_count(self):
		return len(self.get_alignments())

	def get_best_alignments(self, iterate=False):
		return filter(lambda aln: aln.is_best, self.get_alignments())

	def get_best_alignments_count(self):
		return len(self.get_best_alignments())

class Genome(TargetSeq):
	@classmethod
	def from_header(cls, header, length=-1, detailed=False):
		data = db.parse_genome_header(header)
		accession = data.get('ref', 'no-accession')
		gi = int(data['gi'])
		return cls(accession, gi, length, data['description'], detailed=detailed)

	def set_taxid(self, data_access):
		tax = data_access.get_taxids((self.gi,), 'nucl', tuple)
		if tax:
			self.tax_id = tax[0]



class CDS(TargetSeq):
	__slots__ = ('cds', 'nucl_gi',  'prot_gi', 'gene', 'product', 'protein_id')

	@classmethod
	def from_header(cls, cds_header, length, detailed=False):
		data = db.parse_cds_header(cds_header)
		cds = '%s_%s' % (data['gb'], data['cds'])
		tax_id = int(data['taxon'])
		nucl_gi = int(data['nucl_gi'])
		prot_gi = int(data.get('prot_gi', -1))
		protein_id = data.get('protein_id', '')
		gene = data.get('gene', '')
		product = data.get('product', '')
		return cls(cds, tax_id, nucl_gi, length, prot_gi, protein_id,
				   gene, product, detailed=detailed)

	def __init__(self, cds, tax, nucl_gi, length, prot_gi=-1,
				 protein_id='', gene='', product='', detailed=False):
		super(CDS, self).__init__(cds, nucl_gi, length)
		self.cds = cds
		self.tax_id = int(tax)
		self.nucl_gi = int(nucl_gi)
		self.length = int(length)
		self.prot_gi = prot_gi
		self.protein_id = protein_id
		self.gene = gene
		self.product = product
		if detailed:
			self.alignment = DetailedAlignment(self.length)
		else:
			self.alignment = SimpleAlignment(self.length)

	def id(self):
		return self.cds

	def get_key_attributes(self, assigned=False):
		attribs = ['cds', 'tax_id', 'prot_gi', 'length',
			       'protein_id', 'gene', 'product']
		if assigned:
			attribs.append('assigned_total_coverage')
			attribs.append('assigned_coverage_fold')
		else:
			attribs.append('total_coverage')
			attribs.append('coverage_fold')
		attribs.append('get_alignments_count')
		attribs.append('get_best_alignments_count')
		return attribs

	def get_protein_id(self):
		self._return_value('protein_id', None)

	def get_product(self):
		self._return_value('product', None)


class Alignment(object):

	def __init__(self, read_id, target_id, qstart, qend, tstart, tend):
		self.read_id = read_id
		self.target_id = target_id
		self.qstart = qstart
		self.qend = qend
		self.tstart = tstart
		self.tend = tend
		self.read_assigned = False
		self.is_chosen = False
		self.is_best = False

	def is_best_aln(object):
		if 'is_best' not in dir(self):
			raise ValueNotSetError('Value not set at this point of the execution.')
		return self.is_best

	def is_read_assigned(self):
		return self.read_assigned

	def is_chosen(self):
		return self.is_chosen

	def get_score(self):
		raise NotImplementedError("Score not defined for abstract Alignment")


class BlastAlignment(Alignment):
	__slots__ = ['read_id', 'target_id', 'qstart', 'qend',
				 'tstart', 'tend', 'e_value', 'bitscore', 'is_best']

	def __init__(self, read_id, target_id, qstart, qend, tstart, tend, e_value, bitscore):
		super(BlastAlignment, self).__init__(read_id, target_id, qstart,
											 qend, tstart, tend)
		self.e_value = float(e_value)
		self.bitscore = float(bitscore)

	def __cmp__(self, other):
		if not isinstance(other, type(self)):
			raise Exception('Cannot compare two object of different type.')
		return self.bitscore - other.bitscore

	def get_score(self):
		return self.bitscore


class SamAlignment(Alignment):
	__slots__ = ['read_id', 'target_id', 'qstart', 'qend',
				 'tstart', 'tend', 'mapq', 'is_best']
	def __init__(self, read_id, target_id, qstart, qend, tstart, tend, mapq):
		super(SamAlignment, self).__init__(read_id, target_id, qstart,
										   qend, tstart, tend)
		self.mapq = mapq

	def __cmp__(self, other):
		if not isinstance(other, type(self)):
			raise Exception('Cannot compare two object of different type.')
		return self.mapq - other.mapq

	def get_score(self):
		return self.mapq

def create_sequence(header, length, db_type, detailed=False):
	if db_type == 'cds':
		seq = CDS.from_header(header, length, detailed)
	elif db_type == 'genome':
		seq = Genome.from_header(header, length, detailed)
	elif db_type in ('nt', 'nr'):
		seq = TargetSeq.from_header(header, length, detailed)
	return seq

def get_seq_id(header, db_type):
	seq = create_sequence(header, -1, db_type, False)
	return seq.get_id()
