from __future__ import absolute_import
from argparse import ArgumentParser
from datetime import date

import metagenomix
import metagenomix.tax as tax
import metagenomix.utils as utils
import metagenomix.io.datainput as din
import metagenomix.sequence.consensus as consensus_tools


def vcf_header(filters, infos, db_file):
	today = date.today()
	header_seqs = []
	header_seqs.append('##fileformat=VCFv4.2')
	header_seqs.append('##filedate=%s' % today.strftime('%Y%m%d'))
	header_seqs.append('##source=metagenomix_v%s' % metagenomix.__version__)
	header_seqs.append('##reference=%s' % db_file)
	for info in infos:
		header_seqs.append(info.get_vcf_output())
	for filt in filters:
		header_seqs.append(filt.get_vcf_output())
	return '\n'.join(header_seqs)



def get_consensus_parser():
	parser = ArgumentParser(usage="Generate a VCF file and consensus sequence from reference genome and alignment.")
	parser.add_argument(metavar="<REF-DB>", dest="ref_db",
		help="Fasta file with reference sequences")
	parser.add_argument(metavar="<DB-TYPE>", dest="db_type", choices=metagenomix.supported_db_types,
		help="Database type against which reads were aligned")
	parser.add_argument(metavar="<ALN-FILE>", dest="aln_file",
		help="Alignment file")
	parser.add_argument(metavar="<OUTPUT-VFC-FILE>", dest="output_vcf",
		help="File to which VCF will be stored")
	parser.add_argument("--seq", dest="seq_file",
		help="Specify file to which the consensus sequence will be written")
	return parser

def generate_consensus():
	parser = get_consensus_parser()
	args = parser.parse_args()

	if args.seq_file is not None:
		generate_seq, seq = True, args.seq_file
	else:
		generate_seq, seq = False, None

	with utils.timeit('Loading tax tree'):
		tt = tax.TaxTree()

	with utils.timeit('Parsing aln file %s <%s>' % (args.aln_file, utils.get_appropriate_file_size(args.aln_file))):
		read_alns, target_seqs = din.parse(args.aln_file, args.db_type, detailed=True)

	with utils.timeit('Retrieving tax data'):
		utils.retrieve_tax_ids(target_seqs, args.db_type)

	with utils.timeit('Creating consensuses'):
		with open(args.output_vcf, 'w') as fout:
			filters = consensus_tools.get_filters()
			infos = consensus_tools.get_infos()
			fout.write(vcf_header(filters, infos, args.ref_db))
			fout.write('\n')
			fout.write('#CHROM\tPOS\tID\tREF\tALT\tFILTER\tINFO\n')
			for ts in target_seqs.itervalues():
				snps = 0
				for loc, cov, e, fres, infos in consensus_tools.itervariants(ts.alignment, args.output_vcf, generate_seq, seq):
					snps += 1
					filter_output = 'PASS' if fres == [] else ';'.join(fres)
					info_output = ';'.join('%s=%s' % (k, v) for k, v in infos.iteritems())
					fout.write('%s\t%d\t.\t%s\t%s\t%s\t%s\n' % (ts.get_id(), loc, e.target, e.query, filter_output, info_output))
				print '-'*80
				print 'Organism:        ', tt.get_org_name(ts.tax_id)
				print 'Sequence length: ', ts.alignment.length
				print 'Total coverage:  ', ts.alignment.get_coverage()
				print 'Fold:            ', ts.alignment.get_fold()
				print 'Confirmed SNPs:  ', snps
				print