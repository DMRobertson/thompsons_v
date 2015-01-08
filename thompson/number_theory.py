from fractions import gcd

def extended_gcd(a,b):
	"""From `this exposition of the extended gcd algorithm <http://anh.cs.luc.edu/331/notes/xgcd.pdf>`. Computes :math:`d = \gcd(a, b)` and returns a triple :math:`(d, x, y)` where :math:`d = ax + by`."""
	prevx, x = 1, 0; prevy, y = 0, 1
	while b:
		q = a//b
		x, prevx = prevx - q*x, x
		y, prevy = prevy - q*y, y
		a, b = b, a % b
	return a, prevx, prevy

def solve_linear_diophantine(a, b, c):
	r"""Solves the equation :math:`ax + by = c` for integers :math:`x, y \in\mathbb{Z}`.
	
	:rtype: ``None`` if no solution exists; otherwise a triple (base, inc, lcm).
	"""
	d, x, y = extended_gcd(a, b)
	if c % d != 0:
		#Solution exists iff d divides c
		return None
	scale = c // d
	x *= scale
	y *= scale
	assert a*x + b*y == c
	base = (x, y)
	inc = (b // d, a // d)
	lcm = a * b // d
	return base, inc, lcm

def lcm(a, b=None):
	if b is None:
		if len(a) == 0:
			raise ValueError('If computing the LCM of an iterable, the iterable must be nonempty.')
		out = 1
		for number in a:
			out *= lcm(out, number)
		return out
	#otherwise
	return a*b // gcd(a, b)

def solve_linear_congruence(coefficients, values, moduli):
	... #todo