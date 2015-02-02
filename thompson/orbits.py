"""
.. testsetup::
	
	from thompson.orbits import *
	from thompson.examples import *
"""

from collections import namedtuple

from .number_theory import SolutionSet
from .word import format
from .generators import Generators

__all__ = ["ComponentType", "print_component_types", "SolutionSet"]

BaseComponentType = namedtuple("BaseComponentType", "type characteristic")
class ComponentType(BaseComponentType):
	r"""This class is essentially a glorified :mod:`enumeration <py3:enum>`: its instances store a number from 1 to 5 to represent the type of a component.
	
	:ivar characteristic: Either the characteristic :math:`(m, \Gamma)` of this component, or ``None`` if this component has no characteristic. Periodic components have characteristic :math:`(m, \varepsilon)`, where :math:`\varepsilon` is the empty word.
	
	.. seealso:: Section :paperref:`SNF` of the paper.
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
	r"""Prints the classification of the components under *aut* of each word in *words* with respect to *basis*. If *basis* is omitted, it is taken to be the smallest possible expansion which could potentially be a semi-normal basis; see :meth:`~thompson.mixed.MixedAut._seminormal_form_start_point`. If *words* is omited, it is taken to be the same as *basis*.
	
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
