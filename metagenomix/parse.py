class BLASTParser (object):
    ''' Enables BLAST to input format parsing
    '''
    def __init__ (self):
        ''' @param output_format: BLAST output format.
            If none provided, standard output is taken.
            Standard output format is:
            qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        '''
        output_format = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        self._set_up_parser(output_format)

    def _set_up_parser (self, output_format):
        ''' Takes output format specification and sets up
            the parser for the specific format parsing
        '''
        self.fmt_values = {}
        values = output_format.strip().split()
        i = 0
        for value in values:
            self.fmt_values[value] = i
            i += 1

    def parse_file (self, blast_output_fname):
        with open(blast_output_fname) as fin:
            for line in fin:
                line = line.strip()
                if line.startswith('#'):
                    continue
                yield self.parse_line(line)

    def parse_line (self, line):
        aln_data_list   = line.split()
        subject_id      = aln_data_list[self.fmt_values['sseqid']]
        (gi, db_source, nucl_accession) = subject_id.split('|')[1:4]

        score       = aln_data_list[self.fmt_values['bitscore']]
        start       = aln_data_list[self.fmt_values['sstart']]
        stop        = aln_data_list[self.fmt_values['send']]
        query_id    = aln_data_list[self.fmt_values['qseqid']]

        aln_data    = AlignmentData(nucl_accession, db_source, gi,
                                    score, start, stop)
        return (query_id, aln_data)


class AlignmentData (object):

    def __init__ (self, nucleotide_accession, db_source, gi, score, start, stop):
        self.nucleotide_accession   = nucleotide_accession
        self.db_source              = db_source
        self.gi                     = int (gi)
        self.score                  = float (score)
        self.start                  = int(start)
        self.stop                   = int(stop)
        if start > stop:
            self.strand = '+'
            tmp = self.start
            self.start = self.stop
            self.stop = tmp
        else:
            self.strand = '-'

    def __str__(self, *args, **kwargs):
        return "{0},{1},{2},{3},{4},{5},{6}".format (self.nucleotide_accession,
                                                     self.db_source, self.gi,
                                                     self.score, self.start,
                                                     self.stop, self.strand)

def main():
    import sys
    if len(sys.argv) < 2:
        print 'Usage: parse <BLAST-OUTPUT-FILE>'
        sys.exit(-1)

    blast_output = sys.argv[1]
    bp = BLASTParser()
    for read_id, aln_data in bp.parse_file(blast_output):
        print read_id, aln_data

if __name__ == '__main__':
    main()