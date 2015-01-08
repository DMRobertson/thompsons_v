"""
.. testsetup::
	
	from thompson.orbits import *
	from thompson.examples import *
"""

from collections import namedtuple
from numbers import Number

from .number_theory import extended_gcd, solve_linear_diophantine
from .word import format
from .generators import Generators

__all__ = ["ComponentType", "print_component_types", "SolutionSet"]

BaseComponentType = namedtuple("BaseComponentType", "type characteristic")
class ComponentType(BaseComponentType):
	r"""This class is essentially a glorified :mod:`enumeration <py3:enum>`: its instances store a number from 1 to 5 to represent the type of a component.
	
	:ivar characteristic: Either the characteristic :math:`(m, \Gamma)` of this component, or ``None`` if this component has no characteristic. Periodic components have characteristic :math:`(m, \varepsilon)`, where :math:`\varepsilon` is the empty word.
	
	.. seealso:: Section 4.1 of the paper.
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
		"""Describes a complete finite (i.e. periodic, type 2) component. The argument is the period of the the component."""
		return cls(cls._complete_finite, (period, tuple()))
	
	@classmethod
	def semi_infinite(cls, characteristic, backward=False):
		"""Describes a semi-infinite component (types 3, 4). The first argument is the characteristic of the component; use ``None`` if this component doesn't have a characteristic. The second argument should specify where this component is left or right semi-infinite."""
		if characteristic is not None and len(characteristic[1]) == 0:
			raise ValueError('Semi-infinite components cannot have a trivial characteristic multiplier.')
		if characteristic is not None:
			assert (characteristic[0] < 0) == backward, (characteristic, backward)
		type = cls._left_semi_infinite if backward else cls._right_semi_infinite
		return cls(type, characteristic)
	
	@classmethod
	def complete_infinite(cls):
		"""Describes a component which is infinite in both directions (type 1)."""
		return cls(cls._complete_infinite, None)
	
	@classmethod
	def incomplete(cls):
		"""Describes a component which is incomplete finite (type 5)."""
		return cls(cls._incomplete, None)
	
	def __str__(self):
		output = "{} component".format(self._names[self.type])
		if self.type is self._complete_finite:
			output += " of order {}".format(self.characteristic[0])
		if self.type == self._left_semi_infinite or self.type == self._right_semi_infinite:
			if self.characteristic is not None:
				output += " with characteristic ({}, {})".format(
				  self.characteristic[0], format(self.characteristic[1]))
			else:
				output += " with no characteristic"
		return output
	
	def is_type_A(self):
		"""Returns true if this component belongs to an orbit of type A (periodic)."""
		return self.characteristic is not None and len(self.characteristic[1]) == 0
	
	def is_type_B(self):
		"""Returns true if this component belongs to an orbit of type B (has a characteristic)"""
		return self.characteristic is not None and len(self.characteristic[1]) > 0 
	
	def is_type_C(self):
		"""Returns true if this component belongs to an orbit of type C (does not have a characteristic)"""
		return self.characteristic is None
	
	def is_incomplete(self):
		"""Returns True if this component is incomplete, otherwise False."""
		return self.type is self._incomplete

def print_component_types(aut, basis=None, words=None):
	r"""Prints the classification of the components under *aut* of each word in *words* with respect to *basis*. If *basis* is omitted, it is taken to be the smallest possible expansion which could potentially be a semi-normal basis; see :meth:`~thompson.automorphism.Automorphism._seminormal_form_start_point`. If *words* is omited, it is taken to be the same as *basis*.
	
		>>> print_component_types(arity_three_order_inf)
		x1 a1: Left semi-infinite component with characteristic (-1, a1)
		x1 a2: Bi-infinite component
		x1 a3 a1: Bi-infinite component
		x1 a3 a2: Bi-infinite component
		x1 a3 a3: Right semi-infinite component with characteristic (1, a1)
		with respect to the basis [x1 a1, x1 a2, x1 a3 a1, x1 a3 a2, x1 a3 a3]
		
		>>> print_component_types(arity_four)
		x1 a1 a1: Periodic component of order 4
		x1 a1 a2: Periodic component of order 4
		x1 a1 a3: Periodic component of order 4
		x1 a1 a4: Periodic component of order 4
		x1 a2: Bi-infinite component
		x1 a3: Right semi-infinite component with characteristic (1, a3)
		x1 a4: Left semi-infinite component with characteristic (-1, a4)
		with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a1 a3, x1 a1 a4, x1 a2, x1 a3, x1 a4]
	
		.. seealso:: The :meth:`~thompson.automorphism.Automorphism.orbit_type` method.
	"""
	if basis is None:
		basis = aut._seminormal_form_start_point()
	if words is None:
		words = basis
	elif not isinstance(words, Generators):
		words = Generators(aut.domain.signature, words)
	for w in words:
		print("{}: {}".format(
		  w, aut.orbit_type(w, basis)[0]))
	print('with respect to the basis', basis)

BaseSolutionSet = namedtuple('BaseSolutionSet', 'base increment')
class SolutionSet(BaseSolutionSet):
	r"""Create solution sets as follows:
	
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
