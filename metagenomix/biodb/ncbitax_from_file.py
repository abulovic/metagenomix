import os

def loadGi2Taxid(gi2taxid_dump):
    '''
    Maps GI (gene index) values to taxonomical IDs.
    The file containing these values should be obtained
    from the NCBI taxonomy FTP site:
    ftp://ftp.ncbi.nih.gov/pub/taxonomy

    :param gi2taxid_dump path to gi_taxid_[nucl/prot] file.
    Can be obtained from NCBI taxonomy FTP site.
    :rtype dict(key=gi:int, value=taxid:int)
    '''
    if not os.path.isfile(gi2taxid_dump):
        raise ValueError('''Path you supplied to the gi2taxid\
             dump file seems to be invalid.''')
    gi2taxid = {}

    gi2taxid_file = open(gi2taxid_dump, 'r')
    while (True):
        line = gi2taxid_file.readline()
        line = line[0:-1]
        if not line:
            break
        try:
            (gi, taxid) = line.split()
        except ValueError:
            raise ValueError('''Cannot unpack splitted string into\
                two values (gi, taxid) for line %s of input file.'''
                % line)
        gi2taxid[int(gi)] = int(taxid)
    gi2taxid_file.close()

    return gi2taxid

def loadNcbiNames(names_dump):
    '''
    Loads scientific names from NCBI names taxonomy dump.

    :param names_dump path to names dump file
    :rtype dict(key=taxid:int, value=organims_name:str)
    '''
    pass

def loadNcbiRanks(nodes_dump):
    '''
    Loads taxonomy rank from NCBI nodes dump

    :param nodes_dump path to nodes dump file
    :rtype dict(key=taxid:int, value=taxonomy_rank:str)
    '''
    pass


if __name__ == '__main__':
    
    import guppy.heapy.RM
    file_location = raw_input('Enter file location:')
    gi2taxid = loadGi2Taxid(file_location)
    raw_input('Enter to finish off programme:')