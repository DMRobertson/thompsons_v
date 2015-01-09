r"""Functions implementing various number-theoretic algorithms.
 
.. testsetup::
	
	from thompson.examples import *
	from thompson.number_theory import *
"""

from collections import namedtuple
from fractions import gcd
from numbers import Number

__all__ = ['lcm', 'gcd', 'extended_gcd', 'solve_linear_diophantine', 'solve_linear_congruence']

BaseSolutionSet = namedtuple('BaseSolutionSet', 'base increment')
class SolutionSet(BaseSolutionSet):
	r"""Solution sets to the orbit sharing test :math:`u\phi^k = v`. These are either
	
	- empty (no such :math:`k`) exists;
	- a singleton; or
	- a linear sequence a + b*n.
	
	Create solution sets as follows:
	
	.. doctest::
		:options: -ELLIPSIS
		
		>>> print(SolutionSet(5, 2))
		{..., 3, 5, 7, 9, 11, 13, ...}
		>>> print(SolutionSet.singleton(4))
		{4}
		>>> print(SolutionSet.empty_set())
		{}
		>>> print(SolutionSet.the_integers())
		{..., -1, 0, 1, 2, 3, 4, ...}
	"""
	__slots__ = ()
	
	def __new__(cls, base, increment):
		if increment is not None:
			increment = abs(increment)
		self = super(SolutionSet, cls).__new__(cls, base, increment)
		return self
	
	@classmethod
	def singleton(cls, value):
		return cls(value, None)
	
	@classmethod
	def empty_set(cls):
		return cls(None, None)
	
	@classmethod
	def the_integers(cls):
		return cls(0, 1)
	
	def is_sequence(self):
		"""Returns True if this set contains more than one distinct element; otherwise returns False.
		
			>>> SolutionSet.empty_set().is_sequence()
			False
			>>> SolutionSet.singleton(2).is_sequence()
			False
			>>> SolutionSet(base=4, increment=3).is_sequence()
			True
			>>> SolutionSet.the_integers().is_sequence()
			True
		"""
		return self.base is not None and self.increment is not None
	
	def is_singleton(self):
		"""Returns True if this contains precisely one element, otherwise False.
		
			>>> SolutionSet.empty_set().is_singleton()
			False
			>>> SolutionSet.singleton(2).is_singleton()
			True
			>>> SolutionSet(base=4, increment=3).is_singleton()
			False
			>>> SolutionSet.the_integers().is_singleton()
			False
		"""
		return self.base is not None and self.increment is None
	
	def is_empty(self):
		"""Returns True if this set is empty, otherwise False.
		
			>>> SolutionSet.empty_set().is_empty()
			True
			>>> SolutionSet.singleton(2).is_empty()
			False
			>>> SolutionSet(base=4, increment=3).is_empty()
			False
			>>> SolutionSet.the_integers().is_empty()
			False
		"""
		return self.base is None
	
	def is_the_integers(self):
		"""Returns True if this set contains every integer; otherwise returns False.
		
			>>> SolutionSet.empty_set().is_the_integers()
			False
			>>> SolutionSet.singleton(2).is_the_integers()
			False
			>>> SolutionSet(base=4, increment=3).is_the_integers()
			False
			>>> SolutionSet.the_integers().is_the_integers()
			True
		"""
		return self.is_sequence() and abs(self.increment) == 1
	
	def __contains__(self, other): #other in self
		"""Returns true if this set contains an *other* number.
		
			>>> 1024 in SolutionSet.empty_set()
			False
			>>> 1024 in SolutionSet.singleton(128)
			False
			>>> 1024 in SolutionSet(0, 256)
			True
			>>> 1024 in SolutionSet.the_integers()
			True
		"""
		if not isinstance(other, Number):
			return NotImplemented
		
		if self.is_empty():
			return False
		if self.is_singleton():
			return other == self.base
		return self.base % self.increment == other % self.increment
	
	def __and__(self, other): #self = self & other
		"""The ``&`` operator (usually used for bitwise and) stands for intersection of sets.
		
		.. doctest::
			:options: -ELLIPSIS
			
			>>> phi = SolutionSet.empty_set()
			>>> Z = SolutionSet.the_integers()
			>>> singleton = SolutionSet.singleton
			>>> print(phi & phi)
			{}
			>>> print(phi & singleton(1))
			{}
			>>> print(phi & SolutionSet(2, 3))
			{}
			>>> print(singleton(1) & singleton(1))
			{1}
			>>> print(singleton(1) & singleton(2))
			{}
			>>> print(singleton(8) & SolutionSet(4, 2))
			{8}
			>>> print(SolutionSet(1, 3) & SolutionSet(2, 3))
			{}
			>>> print(SolutionSet(1, 3) & SolutionSet(1, 2))
			{..., -5, 1, 7, 13, 19, 25, ...}
			>>> print(SolutionSet(1, 18) & SolutionSet(5, 24))
			{}
			>>> print(SolutionSet(1, 18) & SolutionSet(13, 24))
			{..., -107, -35, 37, 109, 181, 253, ...}
			>>> print(SolutionSet(1, 3) & Z)
			{..., -2, 1, 4, 7, 10, 13, ...}
			>>> print(Z & Z)
			{..., -1, 0, 1, 2, 3, 4, ...}
		"""
		if not isinstance(other, SolutionSet):
			return NotImplemented
		
		if self.is_empty() or other.is_empty():
			return SolutionSet.empty_set()
		
		if self.is_the_integers():
			return other
		
		if other.is_the_integers():
			return self
		
		if self.is_singleton():
			if self.base in other:
				return self
			return SolutionSet.empty_set()
		
		if other.is_singleton():
			if other.base in self:
				return other
			return SolutionSet.empty_set()
		
		#Solve s.base + s.inc * x = o.base + o.inc * y for (x,y)
		a = self.increment
		b = -other.increment
		c = other.base - self.base
		#Equation above is equiv to ax + by = c --- a linear Diophantine equation.
		result = solve_linear_diophantine(a, b, c)
		if result is None:
			return SolutionSet.empty_set()
		base, inc, lcm = result
		new_base = self.base + base[0] * self.increment
		assert new_base in self and new_base in other
		return SolutionSet(new_base, lcm)
	
	def __str__(self):
		if self.is_empty():
			return '{}'
		if self.is_singleton():
			return '{{{0}}}'.format(self.base)
		values = (self.base + i * self.increment for i in range(-1, 5))
		values = (str(num) for num in values)
		return "{{..., {0}, ...}}".format(", ".join(values))


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

def solve_linear_congruence(a, b, n):
	"""Solves the congruence :math:`ax \equiv b \pmod n`.
	
	:rtype: a :class:`SolutionSet`.
	"""
	d = gcd(a, n)
	if b % d == 0:
		return SolutionSet.empty_set()
	a /= d
	b /= d
	n /= d
	
	extended_gcd

