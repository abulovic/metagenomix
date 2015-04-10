import math
import subprocess

_interesting_modules = ('overrepresented sequences', 'per sequence quality scores')


def run_preprocessing(namespace):
	''' Runs all the preprocessing steps - quality estimation, trimming, etc.

	Based on the input parameters, performs the following steps:
	1.  Checks if the FastQC report is available.
	1.a If it is, loads the existing data.
	1.b If it is not, creates a FastQC report in the output directory.
	2.  Performs the quality and adapter trimming (if FastQC detected adapters and contaminant file available)
	3.  Estimates the sample evenness based on sequence overrepresentation in the sample.

	Args:
		namespace: Object (from configparser.parse_args) which contains all the command-line options
	'''
	# Preliminary data generation
	generate_report = True
	input_fname = namespace.fastq_file.split('/')[-1]
	input_prefix = input_fname.split('.')[0]

	# Determine if filtering is needed, and if it is, run it.
	if namespace.fastqc_report:
		generate_report = False
		fastqc_fname = '%s/fastqc_data.txt' % namespace.fastqc_report
	if generate_report:
		cmd = 'fastqc -o %s %s --extract' % (namespace.output_dir, namespace.fastq_file)
		subprocess.call(cmd, shell=True)
		fastqc_fname = '%s/%s_fastqc/fastqc_data.txt' % (namespace.output_dir, input_prefix)

	fastqc_modules = dict(iter(parse_fastqc(fastqc_fname)))
	# Let's stipulate something about the sample evenness and 16s abundance
	overrepresented = fastqc_modules['overrepresented sequences']
	per_seq_qual = fastqc_modules['per sequence quality scores']
	evenness, potential_16s = determine_sample_complexity(overrepresented)
	# Temporary solution
	# TODO
	with open('%s/.preprocessing' % namespace.output_dir, 'w') as fout:
		fout.write('evenness %s\npotential_16s %s\n' % (evenness, potential_16s))
	with open('%s/per-sequence-quality-scores.tsv' % namespace.output_dir, 'w') as fout:
		fout.write('quality	count\n')
		fout.write('\n'.join(map(lambda i: '%d	%s' % (i[0], i[1]), per_seq_qual.items())))

	# Filtering using trim-galore
	filtering = True if overrepresented.seq_constructs else False
	cmd = "trim_galore -q 15 --length 20 "
	if filtering:
		if not namespace.contaminants:
			print 'ERROR: Cannot filter adapters, no contaminant file present. Will perform default cutadapt filtering.'
		else:
			contaminants = overrepresented.get_original_adapters(namespace.contaminants)
			adapters = ' '.join(map(lambda seq: '-a %s' % seq, contaminants.values()))
			cmd += adapters
	cmd = '%s %s > %s/%s.fastq' % (cmd, namespace.fastq_file, namespace.output_dir, input_prefix)
	subprocess.call(cmd, shell=True)
	cmd = 'mv %s_trimmed.fq %s' % (input_prefix, namespace.output_dir)
	subprocess.call(cmd, shell=True)
	cmd = 'mv %s_trimming_report.txt %s' % (input_fname, namespace.output_dir)
	subprocess.call(cmd, shell=True)


def determine_sample_complexity(overrepresented):
	'''Guesses sample complexity based on overrepresented sequences.

	It uses the information on number and percentage of overrepresented sequences
	and tries to guess whether these is rRNA content present in the sample and
	the level of the evenness of the sample:

	Args:
		overrepresented: :metagenomix.preprocessing.fastqc.Overrepresented
	Returns:
		tuple: (evenness, potential_16s), values of which can be low, medium and high.
	'''
	evenness = None
	potential_16s = None
	perc_data = overrepresented.dataset_specific.values()

	# The sample is even and there is low probability of rRNA content.
	if len(perc_data) == 0:
		evenness = 'high'
		potential_16s = 'low'
		return evenness, potential_16s

	# Ignore this for now, it is ashamedly silly:
	N = len(perc_data)
	max_perc = max(perc_data)
	min_perc = min(perc_data)

	if N > 10 and N < 20:
		if max_perc > 5.:
			potential_16s = 'high'
		else:
			potential_16s = 'medium'
	elif N > 20:
		potential_16s = 'high'
	else:
		potential_16s = 'low'

	if max_perc > 0.5 and max_perc < 5.:
		evenness = 'medium'
	elif max_perc > 5.:
		evenness = 'low'
	else:
		evenness = 'high'

	return evenness, potential_16s


class Overrepresented(object):
	'''Place-holder for overrepresented sequences module data.'''
	__slots__ = {'seq_constructs', 'dataset_specific'}

	def __init__(self, seq_constructs=set(), dataset_specific={}):
		self.seq_constructs = seq_constructs
		self.dataset_specific = dataset_specific

	def get_original_adapters(self, contaminant_file):
		contaminants = load_contaminants(contaminant_file)
		return dict(map(lambda n: (n, contaminants[n]), self.seq_constructs))


def _parse_overrepresented_sequences(fin):
	sequencing = set()
	no_hit = {}
	for l in fin:
		l = l.strip()
		if l == '>>END_MODULE':
			return Overrepresented(sequencing, no_hit)
		if l.startswith('#'):
			continue
		seq, count, perc, source = l.split('\t')
		if source.lower() == 'no hit':
			no_hit[seq] = float(perc)
		else:
			parenthesis = source.find('(')
			source = source[: parenthesis - 1]
			sequencing.add(source)

def _parse_per_sequence_quality_scores(fin):
	qual_scores = {}
	for l in fin:
		l = l.strip()
		if l == '>>END_MODULE':
			return qual_scores
		if l.startswith('#'):
			continue
		qual, count = l.split('\t')
		qual_scores[int(qual)] = count

def load_contaminants(contaminant_file):
	'''Loads contaminant sequences from a file.

	Args:
		contaminant_file: path to a file containing sequencing contaminants.
	Returns:
		dict(adapter_name, adapter_sequence)
	'''
	c = {}
	with open(contaminant_file) as fin:
		for l in fin:
			l = l.strip()
			if not l or l.startswith('#'):
				continue
			data = l.split('\t')
			name, seq = data[0], data[-1]
			c[name] = seq
	return c

def parse_fastqc(fqc_fname):
	'''Parses the FastQC report file (fastqc_data.txt in the report directory)

	Function is a generator.
	It checks which of the modules are declared "interesting". For each of the
	interesting modules, invokes a function _parse_<MODULE_NAME>, which returns
	an object with the data contained in the module.
	Overrepresented sequences supported for now.

	Args:
		fqc_fname: Path to the fastqc_data.txt file of the FastQC report.
	Returns:
		Yields tuple (module_name, object with module data)
	'''
	with open(fqc_fname ,'r') as fin:
		for l in fin:
			l = l.strip()
			if l.startswith('#'):
				continue
			if l.startswith('>>') and l != '>>END_MODULE':
				module_name, status = l[2:].lower().split('\t')
				if module_name in _interesting_modules:
					yield module_name, globals()['_parse_%s' % ('_'.join(module_name.split(' ')))](fin)

def main():
	import sys
	if len(sys.argv) < 2:
		print 'Usage: python fasqc.py <FASTQC_DATA.TXT> [CONTAMINANTS_FILE]'
		sys.exit(-1)
	fqc_fname = sys.argv[1]
	modules = {}
	for mod_name, res in parse_fastqc(fqc_fname):
		modules[mod_name] = res
	print modules['overrepresented sequences'].dataset_specific
	print modules['overrepresented sequences'].seq_constructs
	if len(sys.argv) == 3:
		cf = sys.argv[2]
		print modules['overrepresented sequences'].get_original_adapters(cf)



if __name__ == '__main__':
	main()