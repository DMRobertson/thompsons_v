"""Tools for classifying the X<A> orbits of elements of V_n,r.

.. testsetup::
	
	from thompson.orbits import *
"""

from collections import namedtuple
# from enum import Enum

from . import word

__all__ = ["OrbitType", "dump_orbit_types"]

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
		- :math:`z belongs to the orbit being described`; 
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
	"""
	if basis is None:
		basis = aut._minimal_expansion()
	if words is None:
		words = basis
	for w in words:
		print("{}: {}".format(
		  w, aut._orbit_type(w, basis)[0]))
	print('with respect to the basis', basis)
