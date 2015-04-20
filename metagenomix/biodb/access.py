from sqlalchemy.engine import create_engine
from sqlalchemy.orm.scoping import scoped_session
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy import Table, MetaData, Column, Integer, String

from metagenomix.biodb.unity import UnityRecord, UnityCDS
from sqlalchemy.sql.expression import select

_metadata = MetaData()

table_gi_taxid_nuc = Table('nucleotide', _metadata,
    Column('gi', Integer, primary_key=True),
    Column('tax_id', Integer),
)
table_gi_taxid_prot = Table('protein', _metadata,
    Column('gi', Integer, primary_key=True),
    Column('tax_id', Integer),
)
table_ncbi_names = Table('ncbi_names', _metadata,
    Column('tax_id', Integer),
    Column('name_class', String(32)),
    Column('name_txt', String(255))
)

class DbQuery(object):

    supported_tables = ['cds', 'mrna', 'rrna', 'misc_rna']

    '''Serves as a database query utility.'''
    def __init__(self, unity_db_url=None, ncbitax_db_url=None):
        if not unity_db_url:
            unity_db_url = "mysql+mysqldb://root:root@localhost/unity"
        self.unity_db_url = unity_db_url
        if not ncbitax_db_url:
            ncbitax_db_url = "postgres://gi_adm:gi_adm@localhost/gi_taxonomy"
            #ncbitax_db_url = "mysql+mysqldb://root:root@localhost/ncbitax"
        self.ncbitax_db_url = ncbitax_db_url
        self._create_sessions()


    def get_record (self, version, table='cds'):
        '''
        Returns the record associated with the given accession.version.

        :param version: string - GenBank/EMBL/DDBJ/RefSeq Accesion.Version
        :returns: UnityRecord - record associated with the given
                  accession.version. None if no record is found
        '''
        if table not in ('cds', 'rrna', 'mrna', 'misc_rna'):
            raise ValueError('Nonexistent table %s. Only cds, rrna, mrna and misc_rna supported.' % table)

        sess = self.unity_session()
        db_query = """

                SELECT id, db, version, nucl_gi, taxon, location,
                    protein_id, locus_tag, product, gene, prot_gi
                FROM %s
                WHERE version=:version;
            """ % table
        try:

            records = sess.execute(db_query,
            {
                'version': version
             })

            record = None

            for r in records:
                if not record:
                    record = UnityRecord(r['version'])
                cds = UnityCDS(dict(r))
                record.add_cds(cds)

            if record:
                record.cds.sort(key=lambda x: x.location_min)

            return record
        finally:
            self.unity_session.remove()


    def get_taxids (self, gis, table_type='nucl', format=dict):
        '''
        Fetches taxonomy ID for each of the GIs.
        @param gis (list) list of integers representing GIs
        @param format (object type) list or dict.
        @return based on format parameter, returns either list of
        tax IDs or a dictionary mapping gis to tax ids. List can
        contain duplicates.
        '''
        if table_type == 'nucl':
            table = table_gi_taxid_nuc
        elif table_type == 'prot':
            table = table_gi_taxid_prot
        else:
            raise ValueError('%s table type not supported' % table_type)


        if not gis:
            return format()

        sess = self.ncbitax_session()
        try:
            s = select([table.c.gi, table.c.tax_id]
                       ).where(table.c.gi.in_(gis))
            records = sess.execute(s)

            if format == dict:
                gi2taxid_dict = {}
                for (gi, taxid) in records:
                    gi2taxid_dict[int(gi)] = int(taxid)

                return gi2taxid_dict

            elif format == list:
                taxid_list = []
                for (gi, taxid) in records:
                    taxid_list.append (int(taxid))

                return taxid_list

            else:
                return None
        finally:
            self.ncbitax_session.remove()

    def get_organism_names(self, taxids, name_class='scientific name'):
        if not taxids:
            return format()
        sess = self.ncbitax_session()
        try:
            s = select([table_ncbi_names.c.tax_id, table_ncbi_names.c.name_txt]
                       ).where(table_ncbi_names.c.tax_id.in_(taxids)).\
                        where(table_ncbi_names.c.name_class.is_(name_class))
            records = sess.execute(s)

            taxid2name = {}
            for (taxid, name) in records:
                taxid2name[taix] = name
            return taxid2name
        finally:
            self.ncbitax_session.remove()




    def get_organism_name (self, taxid, name_class='scientific name'):
        '''
        Fetches organism name for the speficied taxonomy ID.
        @param taxid (int)  taxonomy ID
        @param name_class (str) scientific name, common name, genbank common
        name, authority
        @return organism name (str)
        '''

        sess = self.ncbitax_session()
        try:

            records = sess.execute("""
                SELECT name_txt
                FROM ncbi_names
                WHERE name_class=:nameClass AND tax_id=:taxId;
            """,
            {
                'nameClass': name_class,
                'taxId': taxid
             })

            record = records.first()

            if record:
                return record['name_txt']

            return None

        finally:
            self.ncbitax_session.remove()

    def get_organism_rank (self, query, by_name=False):
        '''
        Fetches organism rank. Query can be done using organism name
        or organism tax ID.
        @param query (str/int) depends on by_name parameter
        @param by_name (boolean) true if query should be done using organism
        name instead of tax ID.
        @return (str) organism taxonomy rank
        '''
        if by_name:
            tax_id = self.get_organism_taxid(query)
        else:
            tax_id = query
        if not tax_id:
            return None

        tax_id = int(tax_id)

        sess = self.ncbitax_session()
        try:

            records = sess.execute("""
                SELECT rank
                FROM ncbi_nodes
                WHERE tax_id=:taxId;
            """,
            {
                'taxId': tax_id
             })

            record = records.first()

            if record:
                return record['rank']

            return None

        finally:
            self.ncbitax_session.remove()


    def get_organism_taxid (self, organism_name, name_class='scientific name'):
        '''
        Fetches organism taxid for the specified organism name.
        @param organism_name (str) organism nam
        @return taxid (int)
        '''
        sess = self.ncbitax_session()
        try:

            records = sess.execute("""
                SELECT tax_id
                FROM ncbi_names
                WHERE name_class=:nameClass AND name_txt=:nameText;
            """,
            {
                'nameClass': name_class,
                'nameText': organism_name
             })

            record = records.first()

            if record:
                return int(record['tax_id'])

            return None

        finally:
            self.ncbitax_session.remove()


    def _create_sessions(self):
        ''' Creates database sessions '''
        unity_engine = create_engine (self.unity_db_url, echo=False,
                                convert_unicode=True, encoding='utf-8',
                                pool_recycle=3600)
        unity_session = scoped_session(sessionmaker(
                        bind=unity_engine, autocommit=False, autoflush=False))

        self.unity_session = unity_session

        ncbitax_engine = create_engine (self.ncbitax_db_url, echo=False,
                                convert_unicode=True, encoding='utf-8',
                                pool_recycle=3600)

        ncbitax_session = scoped_session(sessionmaker(
                        bind=ncbitax_engine, autocommit=False, autoflush=False))

        self.ncbitax_session = ncbitax_session

class WrongTableError(ValueError):
    def __init__(self, value):
        self.message = "Unsupported table query (%s). Supported tables are: %s" % (value, ','.join(DbQuery.supported_tables))



