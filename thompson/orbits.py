"""Tools for classifying the X<A> orbits of elements of V_n,r.

.. testsetup::
	
	from thompson.orbits import *
"""

from collections import namedtuple
from numbers import Number

from . import word

__all__ = ["OrbitType", "dump_orbit_types", "SolutionSet"]

BaseOrbitType = namedtuple("BaseOrbitType","type data")
class OrbitType(BaseOrbitType):
	r"""Let :math:`y` be a word, and let :math:`X` be an expansion of the standard basis :math:`\boldsymbol{x}=\{x_1, \dotsc, x_n\}`; finally let :math:`\phi` be an  :class:`Automorphism`.
	
	We call the set :math:`\{ y \phi^i\}_{i\in \mathbb Z}` the :math:`\phi`-orbit of :math:`y`. In [Hig]_ (section 9), Higman showed that these :math:`\phi`-orbits could only intersect with :math:`X\langle A \rangle` in certain ways. Here, :math:`X\langle A \rangle` is the set of :meth:`simple <thompson.word.Word.is_simple>` words starting with an :math:`x \in X`.
	
	Let us refer to the intersection of the :math:`\phi`-orbit of :math:`y` with :math:`X\langle A \rangle` as just the *orbit* of :math:`y`. These orbits come in five different types:
	
	1. *Complete infinite.* The orbit exists for all :math:`i \in \mathbb Z` and each element of the orbit is different.
	2. *Complete finite*, or *periodic.* The orbit exists for all :math:`i \in \mathbb Z`, but eventually it repeats itself.
	3. *Right semi-infinite.* The forward orbit :math:`\{ y \phi^n\}_{n\in \mathbb N}` exists and does not repeat itself; however the backward orbit :math:`\{ y \phi^{-n}\}_{n\in \mathbb N}` eventually leaves :math:`X\langle A\rangle`.
	4. *Left semi-infinite.* The backward orbit exists and does not repeat itself; the forward orbit eventually leaves :math:`X\langle A\rangle`.
	5. *Incomplete*. Only a finite number of images :math:`y\phi^{-n}, \dotsc, y\phi^{-1}, y, y\phi, \dotsc, y\phi^m` belong to :math:`X\langle A\rangle`.
	
	This class is essentially a glorified :mod:`enumeration <py3:enum>`: its instances store a number from 1 to 5 to represent the type of an orbit. Instances have a *data* attribute, for storing additional object data. This stores:
	
	* For complete finite orbits, the period of the orbit;
	* For semi-infinite orbits, the characteristic :math:`(m, \Gamma)` of the orbit;
	* For complete infinite orbits, a pair :math:`(z, \Delta)` as described in lemma 4.14(C). Briefly, math:`z` is semi-infinite, and :math:`z\Delta` belonging to the orbit being described;
	* For incomplete orbits, the constant ``None``.
	"""
	_complete_infinite = 1
	_complete_finite = 2
	_right_semi_infinite = 3
	_left_semi_infinite = 4
	_incomplete = 5
	_names = {
		1: "Bi-infinite",
		2: "Periodic",
		3: "Right semi-infinite",
		4: "Left semi-infinite",
		5: "Incomplete"
	}
	
	#Like vanilla namedtuples, do not allow any instance attributes.
	__slots__ = ()
	
	@classmethod
	def periodic(cls, period):
		"""Describes a complete finite (i.e. periodic) orbits. The argument is the period of the the orbit."""
		return cls(cls._complete_finite, period)
	
	@classmethod
	def semi_infinite(cls, characteristic, backward):
		"""Describes a semi-infinite orbit. If backward is True, the orbit is left semi-infinite; otherwise the orbit is right semi-infinite.
		
		:raises ValueError: if direction is neither ``'left'`` nor ``'right'``.
		"""
		if backward:
			type = cls._left_semi_infinite
		else:
			type = cls._right_semi_infinite
		
		assert (characteristic[0] < 0) == (backward)
		return cls(type, characteristic)
	
	@classmethod
	def complete_infinite(cls, type_b_data):
		"""Instance for representing orbits which are infinite in both directions. The argument *type_b_data* is a triple :math:`(n, z, \Delta)` where:
		- :math:`z` belongs to the orbit being described; 
		- we had to move :math:`n` steps forward in the orbit to find :math:`z\Delta`;
		- :math:`z` is a semi-infinite basis element.
		"""
		return cls(cls._complete_infinite, type_b_data)
	
	def __str__(self):
		output = "{} orbit".format(self._names[self.type])
		if self.type is self._complete_finite:
			output += " of order {}".format(self.data)
		elif self.type == self._right_semi_infinite or self.type == self._left_semi_infinite:
			output += " with characteristic ({}, {})".format(
			  self.data[0], word.format(self.data[1]))
		elif self.type == self._complete_infinite:
			output += " containing [{}] {}".format(
			  str(self.data[1]), word.format(self.data[2]))
		return output
	
	def is_type(self, letter):
		"""Returns True if this orbit is of the given type; else False.
		
		:raises ValueError: if *letter* is not one of A, B or C."""
		if letter == "A":
			return self.type == self._complete_finite
		if letter == "B":
			return self.type == self._left_semi_infinite or self.type == self._right_semi_infinite
		if letter == "C":
			return self.type == self._complete_infinite
		raise ValueError("Letter argument should be A, B or C.")

OrbitType.incomplete = OrbitType(OrbitType._incomplete, None)
"""Instance for incomplete orbits."""

#todo renameme
def dump_orbit_types(aut, basis=None, words=None):
	r"""Prints the classification of the orbits under *aut* of each word in *words* with respect to *basis*. If *basis* is omitted, it is taken to be the minimal expansion given by :meth:`~thompson.automorphism._minimal_expansion`. If *words* is omited, it is taken to be the same as *basis*. See the docstring for :meth:`~thompson.automorphism._orbit_type`.
	
		>>> from thompson.examples import *
		>>> dump_orbit_types(arity_three_order_inf)
		x1 a1: Left semi-infinite orbit with characteristic (-1, a1)
		x1 a2: Bi-infinite orbit containing [x1 a3 a3] a2
		x1 a3 a1: Bi-infinite orbit containing [x1 a3 a3] a3
		x1 a3 a2: Bi-infinite orbit containing [x1 a3 a3] a2
		x1 a3 a3: Right semi-infinite orbit with characteristic (1, a1)
		with respect to the basis [x1 a1, x1 a2, x1 a3 a1, x1 a3 a2, x1 a3 a3]
		
		>>> dump_orbit_types(arity_four)
		x1 a1 a1: Periodic orbit of order 4
		x1 a1 a2: Periodic orbit of order 4
		x1 a1 a3: Periodic orbit of order 4
		x1 a1 a4: Periodic orbit of order 4
		x1 a2: Bi-infinite orbit containing [x1 a3] a2
		x1 a3: Right semi-infinite orbit with characteristic (1, a3)
		x1 a4: Left semi-infinite orbit with characteristic (-1, a4)
		with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a1 a3, x1 a1 a4, x1 a2, x1 a3, x1 a4]
	"""
	if basis is None:
		basis = aut._minimal_expansion()
	if words is None:
		words = basis
	for w in words:
		print("{}: {}".format(
		  w, aut._orbit_type(w, basis)[0]))
	print('with respect to the basis', basis)

BaseSolutionSet = namedtuple('BaseSolutionSet', 'base increment')
class SolutionSet(BaseSolutionSet):
	r"""Solutions to the equation :math:`u\psi^m = v` come in specific instances. If there are no solutions, the solution set is :math:`\emptyset`. If :math:`u = v`, then the solution set is :math:`\mathbb{Z}`. Otherwise :math:`u` and :math:`v` share an orbit. If this orbit is periodic, then the solution set is :math:`m + p\mathbb{Z}`, where :math:`p` is the period of the orbit. Otherwise the solution set is a single integer :math:`m`.
	
	Internally we represent this as a pair *(base, increment)* of python objects. The first element *base* is a solution :math:`m` if it exists; otherwise ``None``. The second element *increment* is the increment *p* between solutions (which occurs for a periodic orbit only).
	
	Create solution sets as follows:
	
	.. doctest::
		:options: -ELLIPSIS
		
		>>> print(SolutionSet(5, 2))
		{..., 3, 5, 7, 9, 11, 13, ...}
		>>> print(SolutionSet.singleton(4))
		{4}
		>>> print(SolutionSet.empty_set)
		{}
		>>> print(SolutionSet.all_integers)
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
	
	def is_singleton(self):
		"""Returns True if this contains precisely one element, otherwise False.
		
			>>> SolutionSet.empty_set.is_singleton()
			False
			>>> SolutionSet.singleton(2).is_singleton()
			True
			>>> SolutionSet(base=4, increment=3).is_singleton()
			False
			>>> SolutionSet.all_integers.is_singleton()
			False
		"""
		return self.base is not None and self.increment is None
	
	def is_empty(self):
		"""Returns True if this set is empty, otherwise False.
		
			>>> SolutionSet.empty_set.is_empty()
			True
			>>> SolutionSet.singleton(2).is_empty()
			False
			>>> SolutionSet(base=4, increment=3).is_empty()
			False
			>>> SolutionSet.all_integers.is_empty()
			False
		"""
		return self.base is None
	
	def is_sequence(self):
		"""Returns True if this set contains more than one distinct element; otherwise returns False.
		
			>>> SolutionSet.empty_set.is_sequence()
			False
			>>> SolutionSet.singleton(2).is_sequence()
			False
			>>> SolutionSet(base=4, increment=3).is_sequence()
			True
			>>> SolutionSet.all_integers.is_sequence()
			True
		"""
		return self.base is not None and self.increment is not None
	
	def is_the_integers(self):
		"""Returns True if this set contains every integer; otherwise returns False.
		
			>>> SolutionSet.empty_set.is_the_integers()
			False
			>>> SolutionSet.singleton(2).is_the_integers()
			False
			>>> SolutionSet(base=4, increment=3).is_the_integers()
			False
			>>> SolutionSet.all_integers.is_the_integers()
			True
		"""
		return self.is_sequence() and abs(self.increment) == 1
	
	def __contains__(self, other): #other in self
		"""Returns true if this set contains an *other* number.
		
			>>> 1024 in SolutionSet.empty_set
			False
			>>> 1024 in SolutionSet.singleton(128)
			False
			>>> 1024 in SolutionSet(0, 256)
			True
			>>> 1024 in SolutionSet.all_integers
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
			
			>>> phi = SolutionSet.empty_set
			>>> Z = SolutionSet.all_integers
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
			return SolutionSet.empty_set
		
		if self.is_the_integers():
			return other
		
		if other.is_the_integers():
			return self
		
		if self.is_singleton():
			if self.base in other:
				return self
			return SolutionSet.empty_set
		
		if other.is_singleton():
			if other.base in self:
				return other
			return SolutionSet.empty_set
		
		#Solve s.base + s.inc * x = o.base + o.inc * y for (x,y)
		a = self.increment
		b = -other.increment
		c = other.base - self.base
		#Equation above is equiv to ax + by = c --- a linear Diophantine equation.
		result = solve_linear_diophantine(a, b, c)
		if result is None:
			return SolutionSet.empty_set
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

SolutionSet.all_integers = SolutionSet(base=0, increment=1)
SolutionSet.empty_set = SolutionSet(base=None, increment=None)

def extended_gcd(a,b):
	"""From `this exposition of the extended gcd algorithm <http://anh.cs.luc.edu/331/notes/xgcd.pdf>`. Computes :math:`d = \gcd(a, b)` and returns a triple :math:`(d, x, y)` where :math:`d = ax + by`.
	"""
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

