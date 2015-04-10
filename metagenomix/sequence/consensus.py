'''
A module for generating consensus sequence from the detailed alignment.
Offers functionality to output consensus sequences and VCF files.
'''

import numpy as np
from Bio import SeqIO

from metagenomix import __version__
import metagenomix.utils as utils

_vcf_type_names = {
	int: "Integer",
	float: "Float",
	bool: "Flag",
	str: "String"
}

'''
Filter to evaluate the quality of a variant
Args:
	name: short name of a filter
	description: A long description of what the filter does
	pass_func: A function which accepts a distribution of variants,
			   location and returns a boolean value
'''
class Filter(object):
	def __init__(self, name, description, pass_func):
		self.name = name
		self.description = description
		self.pass_func = pass_func

	def evaluate(self, location, coverage, error):
		return self.pass_func(location, coverage, error)

	def get_vcf_output(self):
		return "FILTER=<ID=%s,Description=\"%s\">" % (self.name, self.description)


def c2(coverage, errors, error):
	if error.cnt >= 3:
		return True
	else:
		return False

c2_filter = Filter('c2', 'Coverage less than 2', c2)

def s50(coverage, errors, error):
	td = total_depth(coverage, errors, error)
	if error.cnt < td/2:
		return False
	else:
		return True

s50_filter = Filter('s50', 'Less than 50%% of samples support call', s50)

def get_filters():
	return [c2_filter, s50_filter]


class Info(object):
	def __init__(self, id, number, type, description, evaluate_func):
		if type not in _vcf_type_names.keys():
			raise ValueError("Info type %s not supported" % str(type))
		self.id = id
		self.number = int(number)
		self.type = type
		self. description = description
		self.evaluate_func = evaluate_func

	def get_vcf_output(self):
		return "INFO=<ID=%s,Number=%d,Type=%s,Description=\"%s\">" % (self.id, self.number,
																	  _vcf_type_names[self.type],
																	  self.description)
	def get_entry(self, cnt):
		return '%s=%d' % (self.id, cnt)

	def evaluate(self, coverage, errors, error):
		return self.evaluate_func(coverage, errors, error)


def supporting_data(coverage, errors, error):
	return error.cnt

def total_depth(coverage, errors, error):
	return coverage + sum([e.cnt for e in errors])

ns_info = Info('NS', 1, int, "Number of samples with data", supporting_data)
dp_info = Info('DP', 1, int, "Total depth", total_depth)

def get_infos():
	return [ns_info, dp_info]

'''
Generates a simple counting consensus from a detailed alignment.
Args:
	alignment: metagenomix.sequence.analysis.DetailedAlignment object
	vcf_file: output file for the consensus
	output_seq: if True, the sequence will be written to a fasta file
	output_fasta: if output_seq is True, needs to be passed to the function
				  for outputing a consensus sequence
'''
def itervariants(alignment, vcf_file, output_seq=False, output_fasta=None):
	if output_seq and output_fasta == None:
		raise ValueError("If output_seq is True, output_fasta argument must be passed.")

	filters = [c2_filter, s50_filter]
	infos = [ns_info, dp_info]
	cov = alignment.coverage

	step = max(1, alignment.length/100)
	valid_mutations = 0
	utils.progressbar(0, start=True)
	for i in sorted(alignment.errors.keys()):
		if i % step == 0:
			utils.progressbar(i/step)
		for e in alignment.errors[i]:
			filt_res = []
			for f in filters:
				result = f.evaluate(cov[i], alignment.errors[i], e)
				if not result:
					filt_res.append(f.name)
			info_res = {}
			for info in infos:
				info_res[info.id] = info.evaluate(cov[i], alignment.errors[i], e)

			if len(filt_res) == 0:
				valid_mutations += 1
				yield i, cov[i], e, filt_res, info_res
	utils.progressbar(100, end=True)
