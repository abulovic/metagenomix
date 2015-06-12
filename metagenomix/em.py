import math
from collections import defaultdict


class Organism(object):
	def __init__(self, tax, unique=0, ambiguous=0):
		self.tax = tax
		self.unique = unique
		self.ambiguous = ambiguous
		self.deltas = {}

	def set_pi(self, pi):
		self.pi = pi

	def set_theta(self, theta):
		self.theta = theta

	def set_mapping_qualities(self, read2quality):
		self.read2quality = {}
		max_qual = max(read2quality.values())
		total_qual = sum([math.exp(qual/max_qual) for qual in read2quality.itervalues()])
		for read, qual in read2quality.iteritems():
			self.read2quality[read] = math.exp(qual/max_qual) / total_qual


# Assigns initial pi and theta estimates to organisms based
# on the number of uniquely and ambiguously mapped reads to
# each of the organisms
def assign_initial_estimates(organisms):
	total_unique = sum([org.unique for org in organisms.itervalues()])
	total_ambiguous = sum([org.ambiguous for org in organisms.itervalues()])

	for org in organisms.itervalues():
		org.set_pi(1./len(organisms))
		#org.set_pi(float(org.unique) / total_unique)
		if total_ambiguous == 0:
			org.set_theta(0.)
		else:
			org.set_theta(float(org.ambiguous) / total_ambiguous)



def exp_max(read_alns, target_seqs, tax_tree):
	present_taxa = set([ts.tax_id for ts in target_seqs.itervalues()])

	err_limit = 1e-20 * math.pow(10, len(present_taxa)/1000)
	all_reads = set(read_alns.keys())
	unique = defaultdict(int)
	ambiguous = defaultdict(int)
	ys = {}

	# we first calculate the hyperparameters (as and bs) by
	# counting the number of uniquely and ambiguously mapped reads
	# to each of the genomes with at least one alignment
	organisms = {}
	for read, alns in read_alns.iteritems():
		if len(alns) == 1:
			ys[read] = 1
			tax = target_seqs[alns[0].target_id].tax_id
			unique[tax] += 1
		else:
			ys[read] = 0
			for tax in [target_seqs[aln.target_id].tax_id for aln in alns]:
				ambiguous[tax] += 1

	# we need to separate alignments to certain taxa:
	tax2alns = defaultdict(dict)
	for read, alns in read_alns.iteritems():
		for aln in alns:
			tax = target_seqs[aln.target_id].tax_id
			tax2alns[tax][read] = aln.get_score()							# TODO! (WILL WORK ONLY FOR BLAST)

	# now we create our data containers
	for tax in present_taxa:
		organisms[tax] = Organism(tax, unique[tax], ambiguous[tax])
		org = organisms[tax]
		org.set_mapping_qualities(tax2alns[tax])

	assign_initial_estimates(organisms)

	#for org in organisms.itervalues():
	#	print tax_tree.get_org_name(org.tax)
	#	print '[UNIQ, AMBIG] = {0}, {1}'.format(org.unique, org.ambiguous)
	#	print '[pi, theta] = {0}, {1}'.format(org.pi, org.theta)

	all_as = sum(unique.values())
	all_bs = sum(ambiguous.values())

	_try = 1
	pi_diff = 0.
	while True:
		print 'TRY', _try
		old_pi_diff = pi_diff
		# we must first calculate denominator sums for each of the reads
		delta_sums = defaultdict(float)
		N = 0.
		for read in all_reads:
			for org in organisms.itervalues():
				if read in org.read2quality:
					org_read = org.pi * math.pow(org.theta, (1-ys[read])) * org.read2quality[read]
					delta_sums[read] += org_read
					N += org_read

		for org in organisms.itervalues():
			for read, qual in org.read2quality.iteritems():
				if delta_sums[read] != 0:
					org.deltas[read] = org.pi * math.pow(org.theta, (1-ys[read])) * org.read2quality[read] / delta_sums[read]
				else:
					org.deltas[read] = 0.


		old_pis = [org.pi for org in organisms.values()]
		for org in organisms.itervalues():
			new_pi = (sum(org.deltas.values()) + unique[org.tax]) / (N + all_as)
			org.set_pi(new_pi)

			new_theta = (sum(map(lambda read: (1 - ys[read]) * org.deltas[read], org.deltas.keys())) + ambiguous[org.tax]) / (sum([1-ys[read] for read in org.deltas.keys()]) + all_bs)
			org.set_theta(new_theta)

		new_pis = [org.pi for org in organisms.values()]
		pi_diff = sum(map(lambda x: (x[0]-x[1])**2, zip(new_pis, old_pis)))
		#for org in organisms.values():
		#	if org.pi > err_limit:
		#		print org.tax, tax_tree.get_org_name(org.tax)
		#		print org.pi
		#		print org.theta
		#		print
		if pi_diff < err_limit:
			break
		_try += 1
		print 'DIFF = {} / {}'.format(pi_diff, err_limit)
		print 'D1-D2 = {}'.format(abs(old_pi_diff - pi_diff))

	reported = {}
	total = 0.
	for tax, org in organisms.iteritems():
		if org.tax == -1 or org.pi < err_limit:
			continue
		total += org.pi
		reported[tax] = org

	for tax, org in reported.iteritems():
		org.perc = org.pi / total
		print tax_tree.get_org_name(tax)
		print org.perc
		print
