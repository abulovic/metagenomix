'''
A module for loading the input data, typically consisting of:
a) a single-end fasta/fastq or paired-end two fasta/fastq files
b) an alignment file corresponding to a certain database and input fasta/fastq files

The module offers loading all the data to RAM, as well as creating and filling
up a database with the data for later reuse.
'''

import os
import json
import sqlite3
from itertools import izip
from contextlib import nested
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metagenomix.biodb.data_access import DataAccess
from metagenomix.parse import BLASTParser
from metagenomix.utils import timeit, time_exec, get_user_reply, table_exists

import metagenomix
import metagenomix.utils as utils
import metagenomix.io.sam as sam
import metagenomix.io.blast as blast

_parse_method = {
	'blast': blast.parse_tab_delimited,
	'sam': sam.parse_cds_sam,
	'xml': blast.parse_xml
}

def parse(aln_file, db_type, *args, **kwargs):
	aln_type = utils.get_file_type(aln_file)
	if aln_type not in metagenomix.supported_aln_types:
		raise ValueError("File type %s not supported." % str(aln_type))
	else:
		return _parse_method[aln_type](aln_file, db_type, *args, **kwargs)


_READ_TABLE_NAME_ = 'read'
_ALN_TABLE_NAME_ = 'alignment'

_json_dict = defaultdict(dict)

def _get_lineage(lineage, json_dict):
	current_dict = json_dict
	for node in lineage[:-1]:
		children = current_dict["children"]
		for child in children:
			if child["id"] == node:
				current_dict = child
	return current_dict

def _traverse(root, present_nodes, tree, abundance):
	if root not in tree.nodes:
		return
	lineage = list(iter(tree.get_lineage(root)))
	children = set(tree.child_nodes[root]) & set(present_nodes)
	if not lineage:
		_json_dict["id"] = root
		_json_dict["name"] = tree.nodes[root].organism_name
		_json_dict["children"] = []
		_json_dict["size"] = abundance[root]
	else:
		d = _get_lineage(lineage, _json_dict)
		if not d.has_key("children"):
			d["children"] = []
		ch = d["children"]
		new_dict = {"id": root, "name": tree.nodes[root].organism_name, "size": abundance[root]}
		if not children:
			new_dict["name"] = "%s (%d)" % (tree.nodes[root].organism_name, abundance[root])
		ch.append(new_dict)
	for child in children:
		_traverse(child, present_nodes, tree, abundance)

@time_exec("EXPORT OTU TO CSV")
def export_OTU_to_csv(read2tax, tax_tree, csv_file):
	taxa = set(read2tax.values())
	with open(csv_file, 'w') as fout:
		for tax in taxa:
			if tax in tax_tree.nodes:
				lineage = map(lambda t: tax_tree.nodes[t].organism_name, tax_tree.get_lineage(tax))
			else:
				lineage = ['unassigned']
			num_reads = len(filter(lambda r: read2tax[r] == tax, read2tax))
			fout.write('%s,%d\n' % (('-').join(lineage), num_reads))

@time_exec("EXPORT OTU TO JSON")
def export_OTU_to_json(read2tax, tax_tree, json_file):
	global _json_dict
	_json_dict = defaultdict(dict)
	present_nodes = set()
	taxa = set(read2tax.values()) - set([0, -1])
	abundance = defaultdict(int)
	for tax in taxa:
		if tax not in tax_tree.nodes:
			continue
		present_nodes.update(tax_tree.get_lineage(tax))
		abundance[tax] = len(filter(lambda r: read2tax[r] == tax, read2tax))
	if 0 in abundance:
		del abundance[0]
	if -1 in abundance:
		del abundance[-1]

	for node in abundance.keys():
		lineage = list(iter(tax_tree.get_lineage(node)))[:-1]
		for anc in lineage:
			abundance[anc] += abundance[node]

	_traverse(tax_tree.root, list(present_nodes), tax_tree, abundance)
	json_string = json.dumps(_json_dict)
	with open(json_file, 'w') as fout:
		fout.write(json_string)


def save_reads_to_fasta(db, cursor, read_ids, fasta_file):
	cursor.execute('SELECT * FROM read WHERE read_id IN (%s)' % (','.join(['?'] * len(read_ids))), tuple(read_ids))
	with open(fasta_file, 'w') as fout:
		for read_id, seq in cursor.fetchall():
			SeqIO.write(SeqRecord(seq=Seq(seq), id=read_id, description=''), fout, 'fasta')


def iter_input_records(fname, format):
	assert(os.path.isfile(fname))
	assert(format in ('fasta', 'fastq'))

	with open(fname) as fin:
		records = SeqIO.parse(fin, format)
		for record in records:
			yield record

def _drop_if_existing(db, cursor, table_name):
	cursor.execute("SELECT name FROM sqlite_master")
	table_names = map(lambda t: t[0], cursor.fetchall())
	if unicode(table_name) in table_names:
		cursor.execute('DROP TABLE %s;' % table_name)
		db.commit()

def create_read_table(db, cursor, store_seq=True, paired=False):
	if not store_seq:
		table_sql = '''CREATE TABLE %s(read_id TEXT PRIMARY KEY)''' % _READ_TABLE_NAME_
	else:
		if not paired:
			table_sql = '''CREATE TABLE %s(read_id TEXT PRIMARY KEY, sequence TEXT)''' % _READ_TABLE_NAME_
		else:
			table_sql = '''CREATE TABLE %s(read_id TEXT PRIMARY KEY, sequence1 TEXT, sequence2 TEXT)''' % _READ_TABLE_NAME_

	cursor.execute(table_sql)

def create_aln_table(db, cursor):
	table_sql = '''CREATE TABLE %s_tmp(read_id TEXT, version TEXT, nucl_gi INTEGER, tax_id INTEGER, score REAL, start INTEGER, end INTEGER, strand TEXT)''' % _ALN_TABLE_NAME_
	cursor.execute(table_sql)

def create_OTU_results_table(db, cursor):
	table_sql = 'CREATE TABLE otu_results(read_id TEXT, assigned_tax INTEGER, species TEXT)'
	cursor.execute(table_sql)

def fill_read_table(db, cursor, fasta1, fasta2, format, pair_end, store_seq):
	if pair_end:
		with nested(open(fasta1), open(fasta2)) as (fin1, fin2):
			records1 = SeqIO.parse(fin1, format)
			records2 = SeqIO.parse(fin2, format)
			for rec1, rec2 in izip(records1, records2):
				if store_seq:
					cursor.execute('INSERT INTO %s(read_id, sequence1, sequence2) VALUES (?,?,?)' % _READ_TABLE_NAME_,
									(rec1.id, str(rec1.seq), str(rec2.seq)))
				else:
					cursor.execute('INSERT INTO %s(read_id) VALUES (?)' % _READ_TABLE_NAME_, (rec1.id,))
	else:
		with timeit('STORING SEQUENCES'):
			for rec in iter_input_records(fasta1, format):
				if store_seq:
					cursor.execute('INSERT INTO %s(read_id, sequence) VALUES (?,?)' % _READ_TABLE_NAME_,
								   (rec.id, str(rec.seq)))
				else:
					cursor.execute('INSERT INTO %s(read_id) VALUES (?)' % _READ_TABLE_NAME_, (rec.id,))
			db.commit()
	db.close()


def fill_blast_aln_table(db, cursor, data_access, aln_file):
	gis = set()
	sql = 'INSERT INTO %s_tmp(read_id, version, nucl_gi, tax_id, score, start, end, strand) VALUES (?,?,?,?,?,?,?,?)' % _ALN_TABLE_NAME_

	with timeit('PARSING BLAST'):
		for read_id, aln_data in BLASTParser().parse_file(aln_file):
			gis.add(aln_data.gi)

	with timeit('LOADING TAXIDS'):
		tax_ids = data_access.get_taxids(list(gis), format=dict)

	with timeit('INSERTING INTO DATABASE'):
		for read_id, aln_data in BLASTParser().parse_file(aln_file):
			tax_id = tax_ids.get(aln_data.gi, -1)
			cursor.execute(sql, (read_id,
								 aln_data.nucleotide_accession,
								 aln_data.gi,
								 tax_id,
								 aln_data.score,
								 aln_data.start,
								 aln_data.stop,
								 aln_data.strand))
	with timeit('REORDERING DATA'):
		cursor.execute('CREATE TABLE alignment(read_id TEXT, version TEXT, nucl_gi INTEGER, tax_id INTEGER, score REAL, start INTEGER, end INTEGER, strand TEXT)')
		cursor.execute('INSERT INTO alignment (read_id, version, nucl_gi, tax_id, score, start, end, strand) SELECT read_id, version, nucl_gi, tax_id, score, start, end, strand FROM alignment_tmp ORDER BY read_id')
		cursor.execute('DROP TABLE alignment_tmp')
		db.commit()

	db.close()

@time_exec("FILLING OTU RESULTS TABLE")
def fill_OTU_results_database(db, cursor, read2tax, no_aln_reads, noassign_tax, tax_tree):
	if table_exists(db, cursor, 'otu_results'):
		cursor.execute('DROP TABLE otu_results')
	insert_sql = 'INSERT INTO otu_results(read_id, assigned_tax, species) VALUES (?,?,?)'

	# read_id TEXT, assigned_tax INTEGER, species TEXT
	create_OTU_results_table(db, cursor)
	all_reads = sorted(list(set(read2tax.keys()) | set(no_aln_reads)))
	for read in all_reads:
		if read in read2tax:
			tax = read2tax[read]
			if tax == noassign_tax:
				cursor.execute(insert_sql, (read, tax, 'unassignable'))
			else:
				cursor.execute(insert_sql, (read, tax, tax_tree.nodes[tax].organism_name))
		else:
			cursor.execute(insert_sql, (read, -1, 'no_aln'))
	db.commit()


def load_input_in_database(db_name, input_file, format, store_seq=True):
	db = sqlite3.connect(db_name)
	cursor = db.cursor()
	if table_exists(db, cursor, _READ_TABLE_NAME_):
		reply = get_user_reply('The table %s seems to already exist in the database %s. Do you wish to override?' % (_READ_TABLE_NAME_, db_name),
							   options=('Y', 'N'))
		if reply == 'n':
			print 'Nothing to do here, data already in database.'
			return db, cursor

		else:
			cursor.execute('DROP TABLE %s' % _READ_TABLE_NAME_)
	create_read_table(db, cursor)
	fasta2 = None
	pair_end = False
	store_seq = True
	fill_read_table(db, cursor, input_file, fasta2, format, pair_end, store_seq)
	return db, cursor

def load_blast_aln_in_database(db_name, aln_file, data_access=None):
	db = sqlite3.connect(db_name)
	cursor = db.cursor()
	if table_exists(db, cursor, _ALN_TABLE_NAME_):
		reply = get_user_reply('The table "%s" seems to already exist in the database %s. Do you wish to override?' % (_ALN_TABLE_NAME_, db_name),
							   options=('Y', 'N'))
		if reply == 'n':
			print 'Nothing to do here, data already in database.'
			return db, cursor

		else:
			cursor.execute('DROP TABLE %s' % _ALN_TABLE_NAME_)
			if table_exists(db, cursor, 'alignment_tmp'):
				cursor.execute('DROP TABLE %s_tmp' % _ALN_TABLE_NAME_)

	if data_access is None:
		data_access = DataAccess()
	create_aln_table(db, cursor)
	fill_blast_aln_table(db, cursor, data_access, aln_file)
	return db, cursor
