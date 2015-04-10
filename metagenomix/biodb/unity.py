# -*- coding: utf-8 -*-
'''
Created on May 27, 2013

@author: marin
'''
import sys
from _location import Location

standard_fields =[
    'id',
    'db'
    'version',
    'nucl_gi',
    'taxon',
    'location',
    'protein_id',
    'locus_tag',
    'product',
    'gene',
    'prot_gi',
    'cds'
]

class UnityRecord(object):
    def __init__(self, version):
        self.version = version
        self.cds = []

    def add_cds(self, cds):
        self.cds.append(cds)

class UnityCDS(object):
    def __init__(self, attributes={}):
        self.attributes = attributes
        self.origin = None
        self.record_id = self.version #Added because of compatibility with older code
        if self.location:
            self.location_min = Location.fast_min_str(self.location)
        else:
            self.location_min = sys.maxint

    def __getattr__(self, name):
        if self.attributes.has_key(name):
            return self.attributes[name]
        elif name in standard_fields:
            return None
        else:
            raise AttributeError("'%s' object has no attribute '%s'",
                self.__class__.__name__, name)
