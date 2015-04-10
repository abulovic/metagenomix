from utils.location import Location

class MockRecord(object):
    ''' Attributes:
        accession:  nucleotide accession
        cdss:       list of coding sequences in the record
        sources:    list of source elemtns in the record
        gi:         genome index
        name:       same as accession? 
        version:    contains accession.version (example AB000181.1)
    '''
    def __init__(self, accession, gi, seq_len):
        self.accession = accession
        self.name = accession
        self.gi = gi
        self.seq_len = seq_len
        self.cdss = [MockCds(self.accession, "1..%d" % seq_len)]
    
    def find_cds(self, location, complement=False, tolerance=0):
        results = []
        for cds in self.cdss:
            if cds.matches(location, complement, tolerance):
                results.append(cds)
        return results
    
class MockCds(object):
    ''' db_xref:    GI:genome_index (string)
        gene:       gene name
        id:         database id
        location:   location string (can be parsed with Location.from_location_str)
        locus_tag:  locus tag
        product:    protein product
        record:     record object (which contains this cds)
        record_id:  record id
        additional_db_xrefs:    additional db_xref identifiers (list of objects)
    '''

    def __init__ (self, accession, location):
        self.location = location
        self.record_id = accession
        self.product = None
        self.locus_tag = None
        self.gene = None
        self.db_xref = None
    

    def matches(self, location, complement, tolerance):
        l1 = Location.from_location_str(self.location, tolerance)
        return l1.intersects(Location.from_location(location, complement))

    def __hash__ (self):
        return hash ((self.record_id, self.location))

    def __eq__ (self, other):
        if (other == None): return False
        return (self.record_id, self.location) == (other.record_id, other.location)

    def __str__(self):
        return "Cds (record_id:" + str(self.record_id) + ", location:" + str(self.location) + ")"



class MockDbQuery (object):

    def __init__ (self, cds_fname):
        '''
        @param: cds_fname(str) Path to a regular FASTA file 
        from which mock database access creates mock records
        which cdss that match the whole record.
        '''
        self.records = {}
        cds_fhandle = open (cds_fname, 'r')
        seq_len = 0
        accession = ''
        record_data = {}
        for line in cds_fhandle.readlines():
            line = line.strip()
            if line.startswith('>'):
                (gi,gi_num,db_source,accession) = line[1:-1].split('|')
            elif not line:
                record = MockRecord (accession, gi, seq_len)
                self.records[accession.split('.')[0]] = record
                seq_len = 0
            else:
                seq_len += len(line)

        cds_fhandle.close()

    def get_record (self, accession, db_source='tst'):
        name = accession.split('.')[0]
        return self.records[name]