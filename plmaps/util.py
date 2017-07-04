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
	"Assume n >= 1 is a power of two".
	return n.bit_length() - 1
