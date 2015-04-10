# -*- coding: utf-8 -*-
'''
Created on May 11, 2013

@author: marin
'''
import re
import sys

#Regular expressions for location parsing
_solo_location = r"[<>]?\d+"
_pair_location = r"[<>]?\d+\.\.[<>]?\d+"
_between_location = r"\d+\^\d+"

_within_position = r"\(\d+\.\d+\)"
_re_within_position = re.compile(_within_position)
_within_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" \
                   % (_within_position,_within_position)
assert _re_within_position.match("(3.9)")
assert re.compile(_within_location).match("(3.9)..10")
assert re.compile(_within_location).match("26..(30.33)")
assert re.compile(_within_location).match("(13.19)..(20.28)")

_oneof_position = r"one\-of\(\d+(,\d+)+\)"
_re_oneof_position = re.compile(_oneof_position)
_oneof_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" \
                   % (_oneof_position,_oneof_position)
assert _re_oneof_position.match("one-of(6,9)")
assert re.compile(_oneof_location).match("one-of(6,9)..101")
assert re.compile(_oneof_location).match("one-of(6,9)..one-of(101,104)")
assert re.compile(_oneof_location).match("6..one-of(101,104)")

assert not _re_oneof_position.match("one-of(3)")
assert _re_oneof_position.match("one-of(3,6)")
assert _re_oneof_position.match("one-of(3,6,9)")

_simple_location = r"\d+\.\.\d+"
_re_simple_location = re.compile(r"^%s$" % _simple_location)
_re_simple_compound = re.compile(r"^(join|order|bond)\(%s(,%s)*\)$"
                                 % (_simple_location, _simple_location))
_complex_location = r"([a-zA-z][a-zA-Z0-9_]*(\.[a-zA-Z0-9]+)?\:)?(%s|%s)" \
                    % (_pair_location, _solo_location)
_re_complex_location = re.compile(r"^%s$" % _complex_location)
_possibly_complemented_complex_location = r"(%s|complement\(%s\))" \
                                          % (_complex_location, _complex_location)
_re_complex_compound = re.compile(r"^(join|order|bond)\(%s(,%s)*\)$"
                                 % (_possibly_complemented_complex_location,
                                    _possibly_complemented_complex_location))

assert _re_simple_location.match("104..160")
assert not _re_simple_location.match("68451760..68452073^68452074")
assert not _re_simple_location.match("<104..>160")
assert not _re_simple_location.match("104")
assert not _re_simple_location.match("<1")
assert not _re_simple_location.match(">99999")
assert not _re_simple_location.match("join(104..160,320..390,504..579)")
assert not _re_simple_compound.match("bond(12,63)")
assert _re_simple_compound.match("join(104..160,320..390,504..579)")
assert _re_simple_compound.match("order(1..69,1308..1465)")
assert not _re_simple_compound.match("order(1..69,1308..1465,1524)")
assert not _re_simple_compound.match("join(<1..442,992..1228,1524..>1983)")
assert not _re_simple_compound.match("join(<1..181,254..336,422..497,574..>590)")
assert not _re_simple_compound.match("join(1475..1577,2841..2986,3074..3193,3314..3481,4126..>4215)")
assert not _re_simple_compound.match("test(1..69,1308..1465)")
assert not _re_simple_compound.match("complement(1..69)")
assert not _re_simple_compound.match("(1..69)")

assert _re_complex_location.match("AL121804:41..610")
assert _re_complex_location.match("AL121804.2:41..610")
assert _re_complex_compound.match("join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)")
assert not _re_simple_compound.match("join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)")
assert _re_complex_compound.match("join(complement(69611..69724),139856..140650)")

#Trans-spliced example from NC_016406, note underscore in reference name:
assert _re_complex_location.match("NC_016402.1:6618..6676")
assert _re_complex_location.match("181647..181905")
assert _re_complex_compound.match("join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)")
assert not _re_complex_location.match("join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)")
assert not _re_simple_compound.match("join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)")
assert not _re_complex_location.match("join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)")
assert not _re_simple_location.match("join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)")

_fast_min_pattern = re.compile('([\w]+:)|([\D_]+)')

def _point(line, tolerance=0):
    '''
    Utility function used to calculate the exact point.

    :param line: string - location string representation
    :param tolerance: int - additional value to add/subtract to fuzzy
                      boundaries
    :returns: int - single point representing the specific location
    '''
    if line[0] == '<':
        return max(int(line[1:])-tolerance, 1)
    elif line[0] == '>':
        return int(line[1:])+tolerance
    else:
        return int(line)


def _loc(line, strand, tolerance=0):
    """FeatureLocation from non-compound non-complement location (PRIVATE).

    Simple examples,

    >>> _loc("123..456", 1000, +1)
    FeatureLocation(ExactPosition(122), ExactPosition(456), strand=1)
    >>> _loc("<123..>456", 1000, strand = -1)
    FeatureLocation(BeforePosition(122), AfterPosition(456), strand=-1)

    A more complex location using within positions,

    """
    try:
        s, e = line.split("..")
    except ValueError:
        assert ".." not in line
        if "^" in line:
            raise LoactionParsingException(
                'Site between two bases not supported')
        elif re.match('[^.]\.[^.]', line):
            raise LoactionParsingException(
                'Single residue location not supported')
        else:
            #e.g. "123"
            s = _point(line, tolerance) if line[0] == '<' else _point(line)
            e = _point(line, tolerance) if line[0] == '>' else _point(line)
            location = Location()
            location.start = s
            location.end = e
            location.strand = strand
            return location

    location = Location()
    location.start = _point(s, tolerance)
    location.end = _point(e, tolerance)
    location.strand = strand
    return location

def _intersects(l1,l2):
    '''
    Utility function used to determine if two locations intersect

    :param l1: First location
    :param l2: Second location
    :returns: bool - True if locations intersect, False otherwise
    '''

    if not l1.sublocations and not l2.sublocations:
        #simple case - the magic intersection formula :)
        return not(l1.start > l2.end or l2.start > l1.end)
    elif not l1.sublocations: #l2 has multiple sublocations
        intersect = False
        for sub_l2 in l2.sublocations:
            intersect = intersect or _intersects(l1, sub_l2)
        return intersect
    elif not l2.sublocations: #l1 has multiple sublocations
        intersect = False
        for sub_l1 in l1.sublocations:
            intersect = intersect or _intersects(sub_l1, l2)
        return intersect
    else: #both have mutiple sublocations
        intersect = False
        for sub_l1 in l1.sublocations:
            for sub_l2 in l2.sublocations:
                intersect = intersect or _intersects(sub_l1, sub_l2)
        return intersect

def _intersection_single(l1,l2):
    '''
    Utility function used to determine intersection range of two location.

    .. note::

    The complement/strand information is not used, only start/stop properties.

    :param l1: First location
    :param l2: Second location
    :returns: Location object representing the intersection, None if no
              intersection is found.
    '''

    if not(l1.start > l2.end or l2.start > l1.end): # If they intersect
        intersection = Location()
        intersection.start = max(l1.start, l2.start)
        intersection.end = min(l1.end, l2.end)
        return intersection
    else:
        return None

def _intersections(l1,l2):
    '''
    Utility function used to determine intersections of two locations.

    .. note:

    The strand information is not used, only start/stop properties.

    :param l1: First location
    :param l2: Second location
    :returns: list of Locations - One location object for each intersection.
              Empty list is returend if no intersections are found.
    '''

    intersections = []

    if l1.complement != l2.complement:
        return intersections
    elif not _intersects(l1, l2):
        return intersections

    if not l1.sublocations and not l2.sublocations:
        #simple case - the magic intersection formula :)
        its = _intersection_single(l1, l2)
        if its:
            its.strand = -1 if l1.complement else 1
            intersections.append(its)
    elif not l1.sublocations: #l2 has multiple sublocations
        for sub_l2 in l2.sublocations:
            its = _intersection_single(l1, sub_l2)
            if its:
                its.strand = -1 if l1.complement else 1
                intersections.append(its)
    elif not l2.sublocations: #l1 has multiple sublocations
        for sub_l1 in l1.sublocations:
            its = _intersection_single(sub_l1, l2)
            if its:
                its.strand = -1 if l1.complement else 1
                intersections.append(its)
    else: #both have mutiple sublocations
        for sub_l1 in l1.sublocations:
            for sub_l2 in l2.sublocations:
                its = _intersection_single(sub_l1, sub_l2)
                if its:
                    its.strand = -1 if l1.complement else 1
                    intersections.append(its)

    return intersections

def _contains(l1, l2):
    '''
    Utility function used to determine if the first location contains the
    second one

    .. note:

    The strand information is not used, only start/stop properties.

    :param l1: First location
    :param l2: Second location
    :returns: bool - True if the first location contains the
              second one, False otherwise
    '''
    if not l1.sublocations and not l2.sublocations:
        if l1.ref != l2.ref:
            return False
        #simple case
        return l1.start <= l2.start and l2.end <= l1.end
    elif not l1.sublocations: #l2 has multiple sublocations
        contains = True
        for sub_l2 in l2.sublocations:
            contains = contains and _contains(l1, sub_l2)
        return contains
    elif not l2.sublocations: #l1 has multiple sublocations
        contains = False
        for sub_l1 in l1.sublocations:
            contains = contains or _contains(sub_l1, l2)
        return contains
    else: #both have mutiple sublocations
        contains = True
        for sub_l2 in l2.sublocations:
            sub_contains = False
            for sub_l1 in l1.sublocations:
                sub_contains = sub_contains or _contains(sub_l1, sub_l2)
            contains = contains and sub_contains
        return contains

class LoactionParsingException(Exception):
    '''
    Exception thrown when problems arise with location parsing
    '''

    pass

class Location(object):

    def __init__(self):
        self.start = None
        self.end = None
        self.ref = None
        self.strand = None
        self.sub_strand = None
        self.operator = None
        self.sublocations = []

    @property
    def complement(self):
        strand = self.strand
        if self.sub_strand:
            strand *= self.sub_strand
        return strand==-1

    def intersects(self, location, use_complement=True):
        '''
        Determines if the the locations intersect with itself

        :param location: Location to search for intersections
        :param use_complement: bool True if strand information is used to
               determine intersections, False otherwise.
        :returns: bool True if locations intersect, False otherwise
        '''
        if use_complement and location.complement != self.complement:
            return False
        return _intersects(self, location)

    def contains(self, location, use_complement=True):
        '''
        Determines if the the given locations is contained within this one

        :param location: Location to check for containment
        :param use_complement: bool True if strand information is used to
               determine containment, False otherwise.
        :returns: bool True if the given location is contained within this one,
               False otherwise
        '''
        if use_complement and location.complement != self.complement:
            return False
        return _contains(self, location)

    def overlaps(self, location, use_complement=True):
        '''
        Determines if the the given location overlaps with this one.

        .. note:: The overlapping is calculated using only the main locations
            start and end positions. No sublocations are involved in the
            calculation. Because of that, two location could be marked as
            overlapping even if they do not share common base pairs.
            Eg. join(1..10,100..200) and join(50..60,80..90) would be mark as
            overlapped.

        :param location: Location to check for overlapping
        :param use_complement: bool True if strand information is used to
               determine overlapping, False otherwise.
        :returns: bool True if the given location overlaps this one,
               False otherwise
        '''
        if use_complement and location.complement != self.complement:
            return False
        else:
            return not(self.start > location.end or location.start > self.end)

    def references(self):
        '''
        Returns a list of references to other records used in this location

        :returns: list List of location references
        '''
        refs = []
        if self.ref:
            refs.append(self.ref)
        refs.extend([l.ref for l in self.sublocations if l.ref])
        return refs

    def min(self):
        '''
        Returns the minimum start value for this location.

        :returns: int Minimum start value of this location
        '''
        if not self.sublocations:
            return self.start
        else:
            m = reduce(lambda x, y: min(x,y),
                           map(lambda x: x.min(),
                               self.sublocations),
                           sys.maxint)
            return m

    @classmethod
    def from_location_str(cls, location_str, tolerance=0):

        #clean
        line = ''.join(location_str.split())
        strand = 1

        #defines the features on the other strand
        if line.startswith('complement'):
            line=line[11:-1]
            strand=-1

        #Special case handling of the most common cases for speed
        if _re_simple_location.match(line):
            #e.g. "123..456"
            s, e = line.split('..')
            location = Location()
            location.start = int(s)
            location.end = int(e)
            location.strand = strand
            return location

        if _re_simple_compound.match(line):
            #e.g. join(<123..456,480..>500)
            i = line.find('(')
            location = Location()
            location.operator = line[:i]
            #we can split on the comma because these are simple locations
            for part in line[i+1:-1].split(','):
                sub_loc = _loc(part, strand, tolerance)
                location.sublocations.append(sub_loc)

            location.start = location.sublocations[0].start
            location.end = location.sublocations[-1].end
            return location

        #Handle the general case with more complex regular expressions
        if _re_complex_location.match(line):
            #e.g. "AL121804.2:41..610"
            if ':' in line:
                location_ref, line = line.split(":")
                location = _loc(line, strand, tolerance)
                location.ref = location_ref
            else:
                location = _loc(line, strand, tolerance)
            return location

        if _re_complex_compound.match(line):
            i = line.find('(')
            location = Location()
            location.operator = line[:i]
            #split on the comma because of positions like one-of(1,2,3) are not
            #supported
            for part in line[i+1:-1].split(','):
                if part.startswith('complement('):
                    part = part[11:-1]
                    if strand == -1:
                        raise LoactionParsingException(
                            'Double complement detected')
                    part_strand = -1
                else:
                    part_strand = 1 #strand
                if ':' in part:
                    part_ref, part = part.split(':')
                else:
                    part_ref = None

                part_loc = _loc(part, part_strand, tolerance)
                part_loc.ref = part_ref
                part_loc.strand = part_strand

                location.sublocations.append(part_loc)

            strands = set(sf.strand for sf in location.sublocations)
            if len(strands)==1:
                sub_strand = location.sublocations[0].strand
            else:
                # i.e. mixed strands
                raise LoactionParsingException('Multiple strands detected')
            location.start = location.sublocations[0].start
            location.end = location.sublocations[-1].end
            location.strand = strand
            location.sub_strand = sub_strand
            return location

        if "order" in line and "join" in line:
            #See Bug 3197
            msg = ('Combinations of "join" and "order" within the same '
                   'location (nested operators) are illegal!')
            raise LoactionParsingException(msg)

        raise LoactionParsingException("Couldn't parse location")



    @classmethod
    def from_location(cls, location_tuple=(0,0), complement=False):
        loc = Location()
        loc.strand = -1 if complement else 1
        loc.start = location_tuple[0]
        if len(location_tuple) != 2:
            loc.end = loc.start
        else:
            loc.end = location_tuple[1]
        return loc

    @classmethod
    def fast_min_str(cls, location_str):
        positions = [int(p) for p in _fast_min_pattern.sub(
                     ' ', location_str).split()]
        return reduce(lambda x, y: min(x,y),
                           positions,
                           sys.maxint)


    def find_intersection (self, location):
        '''
        Computes intersection region between the supplied location and itself
        :param location: Location to search for intersections
        :return: None if the supplied locations don't intersect,
                 otherwise returns the location object representing the
                 intersections
        '''

        intersections = _intersections(self, location)
        if intersections:
            if len(intersections) > 1: #multiple intersections
                intersection = Location()
                intersection.sublocations = intersections
                intersection.strand = -1 if self.complement else 1
                intersection.sub_strand = 1
                intersection.operator = 'join'
                return intersection
            else: #single intersection
                intersection = intersections[0]
                intersection.strand = -1 if self.complement else 1
                return intersection
        else:
            return None

    def __str__(self):
        row = ""

        if self.sublocations:
            row = ",".join(str(l) for l in self.sublocations)
        else:
            if self.start == self.end:
                row = str(self.start)
            else:
                row = '%d..%d' % (self.start, self.end)

        if self.operator:
            row = '%s(%s)' % (self.operator, row)

        if self.complement:
            row = 'complement(%s)' % row

        return row

    def length(self):
        '''
        Calculates the number of bases covered by this location.

        .. note:

        If location contains sublocations, all the lengths of sublocations
        will be included.

        :returns: int - the length of covered bases
        '''
        ret = 0
        if self.sublocations:
            for subloc in self.sublocations:
                ret += subloc.length()
        else:
            ret = self.end - self.start + 1
        return ret
