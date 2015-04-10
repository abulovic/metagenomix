import numpy as np
from collections import deque, defaultdict


class AlignmentLocationError(Exception):
	def __init__(self, tstart, tend, length):
		self.tstart = tstart
		self.tend = tend
		self.length = length

	def __str__(self):
		return "Location (%d, %d) invalid for sequence of length %d" % (self.tstart, self.tend, self.length)


class Alignment(object):
	__slots__ = ('length', 'locations', 'coverage', 'fold')

	def __init__(self, length):
		self.length = length
		self.locations = []
		self.fold = 0


class SimpleAlignment(Alignment):
	__slots__ = ('mismatches',)

	def add_location(self, tstart, tend, mismatches=0):
		if (tstart > self.length or tend > self.length):
			raise AlignmentLocationError(tstart, tend, self.length)
		if tstart > tend:
			tstart, tend = tend, tstart
		self.locations.append((tstart, tend))
		self.fold += abs(tstart - tend) + 1

	def _get_joint_regions(self):
		if len(self.locations) <= 1:
			return self.locations
		pre = list(self.locations)
		pre.sort()
		pre = deque(pre)
		regions = list()
		regions.append(pre.popleft())
		while len(pre):
			r1 = regions.pop()
			r2 = pre.popleft()
			if r1[1] > r2[0]:
				regions.append((r1[0], r2[1]))
			else:
				regions.append(r1)
				regions.append(r2)
		return regions

	def get_covered_length(self):
		regions = self._get_joint_regions()
		return sum([abs(reg[0]-reg[1])+1 for reg in regions])

	def get_coverage(self):
		return self.get_covered_length() / float(self.length)

	def get_fold(self):
		return float(self.fold) / self.get_covered_length()


class AlnError(object):
	__slots__ = ('location', 'target', 'query', 'cnt')

	def __init__(self, location, target, query):
		self.location = location
		self.target = target
		self.query = query
		self.cnt = 1

	def increment(self):
		self.cnt += 1

	def get_target(self):
		return self.target

	def get_query(self):
		return self.query

	def __eq__(self, other):
		return type(other) == type(self) and self.location == other.location

	def __str__(self):
		return '%s (%s -> %s)' % (self.get_type(), self.target, self.query)


class Substitution(AlnError):
	def get_type(self):
		return 'S'
	def get_data(self):
		return [self.target, self.query, self.cnt]

class Insertion(AlnError):
	def __init__(self, location, query):
		super(Insertion, self).__init__(location, '-', query)
	def get_type(self):
		return 'I'
	def get_data(self):
		return [self.query, self.cnt]

class Deletion(AlnError):
	def __init__(self, location, target):
		super(Deletion, self).__init__(location, target, '-')
	def get_type(self):
		return 'D'
	def get_data(self):
		return [self.target, self.cnt]


class Errors(object):
	__slots__ = ('errors',)

	def __init__(self):
		self.errors = defaultdict(list)

	def add(self, error):
		loc = error.location
		obj = None
		for e in self.errors[loc]:
			if isinstance(e, type(error)) and e == error:
				obj = e
				break
		if obj is not None:
			obj.increment()
		else:
			self.errors[loc].append(error)

	def __str__(self):
		s = ''
		for location, ers in self.errors.iteritems():
			s += '%d: ' % location
			for e in ers:
				s += '%s, ' % str(e)
			s += '\n'
		return s

	def __iter__(self):
		return self.errors.iteritems()

	def __len__(self):
		return len(self.errors)

def get_error(location, target_ascii, query_ascii):
	t = chr(target_ascii)
	q = chr(query_ascii)
	if t == '-':
		return Insertion(location, q)
	elif q == '-':
		return Deletion(location, t)
	else:
		return Substitution(location, t, q)


class DetailedAlignment(Alignment):
	__slots__ = ('errors',)

	def __init__(self, length):
		super(DetailedAlignment, self).__init__(length)
		self.coverage = np.zeros(length, 'i')
		self.errors = Errors()

	def add_location(self, tstart, tend, query, target):
		if (tstart > self.length or tend > self.length):
			raise AlignmentLocationError(tstart, tend, self.length)
		if tstart > tend:
			tstart, tend = tend, tstart

		ascii_query = np.frombuffer(query, np.byte)
		ascii_target = np.frombuffer(target, np.byte)
		insert_locations = np.where(ascii_target == 45)[0]
		insert_nucl = ascii_query[insert_locations]
		for i, iloc in enumerate(insert_locations):
			self.errors.add(get_error(tstart + iloc - i - 1, ascii_target[iloc], ascii_query[iloc]))

		no_insert = np.where(ascii_target != 45)[0]
		no_insert_target = ascii_target[no_insert]
		no_insert_query = ascii_query[no_insert]
		diff = no_insert_target - no_insert_query
		is_zero = np.where(diff==0)[0] + (tstart - 1)
		self.coverage[is_zero] += 1
		non_zero = np.where(diff != 0)[0]
		for nz_loc in non_zero:
			self.errors.add(get_error(nz_loc+tstart-1, no_insert_target[nz_loc], no_insert_query[nz_loc]))

	def get_covered_length(self):
		return np.count_nonzero(self.coverage)

	def get_coverage(self):
		return self.get_covered_length() / float(self.length)

	def get_fold(self):
		return sum(self.coverage) / float(self.get_covered_length())


if __name__ == '__main__':
	sa = SimpleAlignment(100)
	sa.add_location(100, 81)
	sa.add_location(1, 20)
	sa.add_location(1,100)
	#sa.add_location(1, 200)
	print sa.get_coverage()
	print sa.fold, sa.get_covered_length()
	print sa.get_fold()

	da = DetailedAlignment(100)
	da.add_location(1,5, 'ACTGC', 'ATTGC')
	da.add_location(91,100, 'AAAAAAAAAA', 'AAAAAAAA-A')
	print da.coverage
	print da.errors


