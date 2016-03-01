"""Utility functions that try to make brute-forcing the power conjugacy problem managable for examples with small bounds."""

def search_pattern(sbound, obound):
	"""An optimistic search pattern which tries to delay expensive computations until as late as possible.
	
		>>> list(search_pattern(3, 5))
		[(1, 1), (2, 1), (2, 2), (1, 2), (3, 1), (3, 2), (3, 3), (2, 3), (1, 3), (1, 4), (2, 4), (3, 4), (1, 5), (2, 5), (3, 5)]
	"""
	assert sbound <= obound
	for i in range(1, sbound + 1):
		yield from corner(i)
	for b in range(sbound + 1, obound + 1):
		for a in range(1, sbound + 1):
			yield a, b

def corner(radius):
	"""Yields the integer coordinates of the top right corner of a square of side length *radius* centered at the origin.
	
		>>> list(corner(3))
		# [(3, 1), (3, 2), (3, 3), (2, 3), (1, 3)]
	"""
	for i in range(1, radius):
		yield radius, i
	yield radius, radius
	for i in reversed(range(1, radius)):
		yield i, radius

def mirrored(iterator):
	for a, b in iterator:
		yield a, b
		yield -a, b

class PowerCollection(dict):
	def __init__(self, aut):
		super().__init__(self)
		self[1] = aut
		self[0] = type(aut).identity(aut.signature)
		self[-1] = ~aut
	
	def __getitem__(self, power):
		try:
			return super().__getitem__(power)
		except KeyError:
			base = self[1 if power > 0 else -1]
			ref = power - 1 if power > 0 else power + 1
			new = self[ref] * base
			inherit_relabellers(new, base)
			self[power] = new
			return new

def inherit_relabellers(target, source):
	target.domain_relabeller = source.domain_relabeller
	target.range_relabeller = source.range_relabeller
