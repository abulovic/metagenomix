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

class DbQuery(object):

    supported_tables = ['cds', 'mrna', 'rrna', 'misc_rna']

    '''Serves as a database query utility.'''
    def __init__(self, gi_taxonomy_url=None):
        if not gi_taxonomy_url:
            gi_taxonomy_url = "postgres://gi_adm:gi_adm@localhost/gi_taxonomy"
        self.gi_taxonomy_url = gi_taxonomy_url
        self._create_sessions()

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

        sess = self.gi_tax_session()
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
            self.gi_tax_session.remove()

    def _create_sessions(self):
        ''' Creates database sessions '''
        gi_tax_engine = create_engine(self.gi_taxonomy_url, echo=False,
            convert_unicode=True, encoding='utf-8',
            pool_recycle=3600)
        gi_tax_session = scoped_session(sessionmaker(
            bind=gi_tax_engine, autocommit=False, autoflush=False))
        self.gi_tax_session = gi_tax_session


class WrongTableError(ValueError):
    def __init__(self, value):
        self.message = "Unsupported table query (%s). Supported tables are: %s" % (value, ','.join(DbQuery.supported_tables))



