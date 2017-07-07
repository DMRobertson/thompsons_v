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

def powers_of_two_between(a, b):
	assert a < b
	
	while power <= b:
		yield term
		power *= 2

	