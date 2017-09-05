from fractions import Fraction
from itertools import tee

def pairwise(iterable):
	#Nicked from itertools documentation
	"s -> (s0,s1), (s1,s2), (s2, s3), ..."
	a, b = tee(iterable)
	next(b, None)
	return zip(a, b)

def increasing_sequence(iterable):
	return all(
		prev < next
		for prev, next in pairwise(iterable)
	)

def all_satisfy(iterable, func):
	return all( func(x) for x in iterable )

def grad(x0, x1, y0, y1):
	return (y1 - y0) / (x1 - x0)
	
def int_power_of_two(n):
	#Assume n >= 1
	return n & (n - 1) == 0

def lerp(x, x0, x1, y0, y1):
	return y0 + grad(x0, x1, y0, y1) * (x - x0)

def ends(indexable):
	return indexable[0], indexable[-1]

def fixed_point(x0, x1, y0, y1):
	"""Where does the line segment mapping [d0, d1] -> [r0, r1] cross the line :math:`y = x`?
	
	:returns: ``None`` if this segment has gradient one, or there is no crossing.
	"""
	m = grad(x0, x1, y0, y1)
	if m == 1:
		return None
	crossing = (y0 - m * x0) / (1 - m)
	if x0 <= crossing <= x1:
		return crossing

def ilog2(n):
	"Assume n >= 1 is a power of two. What is the "
	"""
	If n > 1 is not a power of two, n.bit_length() == ceil( log2(n))
	If n is a power of two, n.bit_length() == log2(n) + 1.
	So n.bit_length() = floor(log2(n) + 1)
	"""
	return n.bit_length() - 1

def small_divisors(n):
	"""Yields the positive divisors d of n > 0 which are smaller in magnitude than their complement n/d."""
	i = 1
	while i * i <= n:
		if n % i == 0:
			yield i
		i += 1

def gradient_roots_dyadic(initial):
	steep = initial > 1
	initial = max(initial.numerator, initial.denominator)
	log = ilog2(initial)
	for d in small_divisors(log):
		if steep:
			yield Fraction(2**d, 1)
		else:
			yield Fraction(1, 2**d)

def dump(aut, name):
	from thompson import Automorphism
	if isinstance(aut, Automorphism):
		print(name, aut.tikz_path(), sep="\n")
	else:
		print(name, aut, sep="\n")
		 
def glue(*auts, cls=None):
	#TODO: put this function in another file. Circular imports are messy
	from .plmap import PL2
	from .cplmap import CPLMap
	if cls is None:
		cls = PL2
	D = list( auts[0].domain )
	R = list( auts[0].range  )
	for prev, next in pairwise(auts):
		if not issubclass(cls, CPLMap):
			assert next.domain[0] == prev.domain[-1]
			assert next.range [0] == prev.range [-1]
		D += next.domain[1:]
		R += next.range[1:]
	
	if issubclass(cls, CPLMap):
		cls._uncycle(D)
		cls._uncycle(R)
	return cls(D,R)
