from metagenomix.biodb.access import DbQuery

class DataAccess (object):
    '''
    Provides access to CDS data and NCBItax data.
    Encapsulates data access, since all data
    can be accessed in two ways:
    * using MySql database (DATABASE type access)
    * by loading raw data from files (FILE type access)
    This method must provide all methods that would be
    available via standard database querying.
    '''

    def __init__(self):
        self._db_access = DbQuery()
        self._gi2taxid_cache = {}

    def clear_cache(self):
        '''
        Clears gi2taxid cache used for faster
        data loding since each get_taxid query  may execute
        a MySql query.
        '''
        self._gi2taxid_cache = {}

    def get_record(self, version, table='cds'):
        '''
        Returns the record associated with the given accession.version.

        :param version: string - GenBank/EMBL/DDBJ/RefSeq Accesion.Version
        :returns: UnityRecord - record associated with the given
                  accession.version. None if no record is found
        '''
        return self._db_access.get_record(version)

    def get_taxids (self, gis, table_type='nucl', format=dict):
        '''
        Fetches taxonomy ID for each of the GIs.

        :param gis (list) list of integers representing GIs
        "param format (object type) list, set or dict.
        :rtype based on format parameter, returns either list of
        tax IDs or a dictionary mapping gis to tax ids. List can
        contain duplicates.
        '''
        requested_gis = {}
        gis = set(gis)
        retrieved_gis = set(self._gi2taxid_cache.keys())
        needed_gis = gis - retrieved_gis
        requested_gis = dict((k, self._gi2taxid_cache[k]) for k in (gis & retrieved_gis))
        if len(needed_gis):
            requested_gis.update(self._db_access.get_taxids(needed_gis, table_type, format=dict))
        for (gi, taxid) in requested_gis.items():
            self._gi2taxid_cache[gi] = taxid
        if format != dict:
            return requested_gis.values()
        else:
            return requested_gis

    def get_organism_name (self, taxid, name_class='scientific name'):
        return self._db_access.get_organism_name(taxid, name_class)

    def get_organism_rank (self, query, by_name=False):
        return self._db_access.get_organism_rank(query, by_name)


def main():
    access = DataAccess()
    print access.get_taxids([53794196], table_type='nucl')
    print access.get_taxids([515507319], table_type='prot')

if __name__ == '__main__':
    main()
