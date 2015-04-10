from collections import Counter, defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metagenomix.io.datainput import save_reads_to_fasta, fill_OTU_results_database
from metagenomix.utils import table_exists, time_exec, timeit, progressbar

def assign_OTU_from_scores(aln_data, tax_tree, forbidden_taxa=[]):
	forbidden_taxa = list(forbidden_taxa)
	forbidden_taxa.extend([0, -1])
	subspecies_taxa = map(lambda data: data[2], aln_data)
	species_taxa = filter(lambda t: t not in forbidden_taxa,
						  map(lambda t: tax_tree.get_parent_with_rank(t, 'species'), subspecies_taxa))
	if not species_taxa:
		return -1
	c = Counter(species_taxa)
	tax = c.most_common(1)[0][0]
	return tax

def _load_alignment_rows(cursor, start, stop, filter_last=True):
	cursor.execute('SELECT rowid, * FROM alignment where rowid BETWEEN ? AND ?', (start, stop))
	data = cursor.fetchall()
	unique_read_ids = set(map(lambda data: data[1], data))
	if filter_last:
		unique_read_ids -= set(data[-1][1])
	for read in unique_read_ids:
		yield read, map(lambda data: data[2:], filter(lambda data: data[1] == read, data))

def _separate_species(read2tax, tax_tree):
	total_aln_reads = len(read2tax.keys())
	tax2read = defaultdict(list)
	ok_species = set()
	bad_species = set()

	for read, tax_id in read2tax.iteritems():
		tax2read[tax_id].append(read)

	for tax, reads in tax2read.iteritems():
		perc = 100. * len(reads) / total_aln_reads
		if perc < 0.1:
			bad_species.add(tax)
		else:
			ok_species.add(tax)
	return ok_species, bad_species

def _perform_initial_OTU_assignment(db, cursor, tax_tree):
	index = 0
	step = 10000
	read2tax = {}
	cursor.execute('SELECT MAX(rowid) from alignment')
	num_rows, = cursor.fetchall()[0]

	progressbar(0, start=True)
	while index < num_rows - step:
		progressbar(index * 100 / num_rows)
		#print '%3d%%' % (index * 100 / num_rows)
		total_batch_reads = 0
		for read, alns in _load_alignment_rows(cursor, index, index + step):
			read2tax[read] = assign_OTU_from_scores(alns, tax_tree)
			total_batch_reads += len(alns)
		if total_batch_reads == 0:
			step *= 2
		index += total_batch_reads
	else:
		for read, alns in _load_alignment_rows(cursor, index, index + step):
			read2tax[read] = assign_OTU_from_scores(alns, tax_tree)
	progressbar(100, end=True)
	return  read2tax

def _filter_bad_species(cursor, read2tax, bad_species, tax_tree):
	questionable_reads = []
	for read, tax in read2tax.iteritems():
		if tax in bad_species:
			questionable_reads.append(read)

	with timeit('EXTRACTING BAD SPECIES READ INFO'):
		cursor.execute('SELECT * FROM alignment WHERE read_id in (%s)' % (','.join(['?']*len(questionable_reads))), questionable_reads)
		data = cursor.fetchall()

	read2aln = defaultdict(list)
	for d in data:
		read2aln[d[0]].append(d[1:])

	new_read2tax = {}
	with timeit('ASSIGNING NEW OTU TO BAD SPECIES READS'):
		reads = set(map(lambda d: d[0], data))
		for read in reads:
			alns = read2aln[read]
			# alns = map(lambda d: d[1:], filter(lambda d: d[0] == read, data))
			new_read2tax[read] = assign_OTU_from_scores(alns, tax_tree, bad_species)

	return new_read2tax


@time_exec('OTU ANALYSIS')
def perform_OTU_analysis_on_db(db, cursor, tax_tree, no_aln_reads_file, basic_stats_file, phylum_level_composition_file):
	cursor.execute('SELECT read_id from read')
	read_ids = map(lambda t: t[0], cursor.fetchall())

	no_aln_reads = get_basic_stats_on_db(db, cursor, no_aln_reads_file, basic_stats_file)
	# 1. Decide on species based on the distribution of species
	#    among the alignments
	read2tax = _perform_initial_OTU_assignment(db, cursor, tax_tree)
	# 2. Estimate which species are really there.
	#    Based on sythetic dataset analysis, 0.1% is the
	#    lower limit of representation for us to be able to say
	#    with some level of certainty that it is there
	ok_species, bad_species = _separate_species(read2tax, tax_tree)
	# 3. Remove the 'bad' species from the reported species.
	#    Either assign them to an already existing species or
	#    declare them unassignable
	new_read2tax = _filter_bad_species(cursor, read2tax, bad_species, tax_tree)
	read2tax.update(new_read2tax)

	fill_OTU_results_database(db, cursor, read2tax, no_aln_reads, -1, tax_tree)
	analyze_phylum_level_composition(read2tax, tax_tree, phylum_level_composition_file)
	return read2tax, bad_species

@time_exec("PHYLUM LEVEL COMPOSITION")
def analyze_phylum_level_composition(read2tax, tax_tree, phylum_level_composition_file):
	tax2read = defaultdict(int)
	for read, tax in read2tax.iteritems():
		if tax == -1:
			continue
		tax2read[tax] += 1

	present_phyla = defaultdict(int)
	for tax, num_reads in tax2read.iteritems():
		parent = tax_tree.get_parent_with_rank(tax, 'phylum')
		present_phyla[parent] += num_reads

	with open(phylum_level_composition_file, 'w') as phylum_fout:
		for tax, num in sorted(present_phyla.items(), key=lambda key: key[1]):
			phylum_fout.write('%s, %d\n' % (tax_tree.nodes[tax].organism_name, num))



@time_exec('BASIC ALIGNMENT STATS')
def get_basic_stats_on_db(db, cursor, no_aln_reads_file, basic_stats_file):
	assert table_exists(db, cursor, 'read')
	assert table_exists(db, cursor, 'alignment')

	cursor.execute('SELECT read_id FROM read')
	reads = map(lambda t: t[0], cursor.fetchall())

	cursor.execute('SELECT read_id, COUNT(*) FROM alignment GROUP BY read_id')
	aln_reads = dict(cursor.fetchall())

	# Let's save the no-alignment reads to fasta if filename provided
	no_aln_reads = set(reads) - set(aln_reads.keys())
	if no_aln_reads_file is not None:
		save_reads_to_fasta(db, cursor, no_aln_reads, no_aln_reads_file)

	avg_num_alns = sum(aln_reads.values()) / float(len(aln_reads.values()))
	c = Counter(aln_reads.values())
	with open(basic_stats_file, 'w') as stats_fout:
		stats_fout.write('BASIC ALIGNMENT STATISTICS\n\n')
		stats_fout.write('Total number of reads: %d\n' % len(reads))
		stats_fout.write('Reads with no alignments: %d (%.3f)\n' % (len(no_aln_reads), float(len(no_aln_reads)) / len(reads)))
		stats_fout.write('Average number of alignments per read: %.3f\n\n' % avg_num_alns)
		stats_fout.write('# of alns, # of reads with this much alignments\n')

		for (num_alns, num_reads) in c.iteritems():
			stats_fout.write('%d, %d\n' % (num_alns, num_reads))

	return no_aln_reads
