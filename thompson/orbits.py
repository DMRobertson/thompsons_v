"""Tools for classifying the X<A> orbits of elements of V_n,r."""

from collections import namedtuple
# from enum import Enum

__all__ = ["OrbitType", "dump_orbit_types"]

BaseOrbitType = namedtuple('BaseOrbitType', 'type power multiplier')
class OrbitType(BaseOrbitType):
	r"""Let :math:`y` be a word, and let :math:`X` be an expansion of the standard basis :math:`\boldsymbol{x}=\{x_1, \dotsc, x_n\}`; finally let :math:`\phi` be an  :class:`Automorphism`.
	
	We call the set :math:`\{ y \phi^i\}_{i\in \mathbb Z}` the :math:`\phi`-orbit of :math:`y`. In [Hig]_ (section 9), Higman showed that these :math:`\phi`-orbits could only intersect with :math:`X\langle A \rangle` in certain ways. Here, :math:`X\langle A \rangle` is the set of :meth:`simple <thompson.word.Word.is_simple>` words starting with an :math:`x \in X`.
	
	Let us refer to the intersection of the :math:`\phi`-orbit of :math:`y` with :math:`X\langle A \rangle` as just the *orbit* of :math:`y`. This class is two things:
	
	- an :mod:`enumeration <py3:enum>` which labels Higman's five orbit types; see section 4.1 of the paper.
		1. *Complete infinite.* The orbit exists for all :math:`i \in \mathbb Z` and each element of the orbit is different.
		2. *Complete finite.* The orbit exists for all :math:`i \in \mathbb Z`, but eventually it repeats itself.
		3. *Right semi-infinite.* The forward orbit :math:`\{ y \phi^n\}_{n\in \mathbb N}` exists and does not repeat itself; however no part of the backward orbit :math:`\{ y \phi^{-n}\}_{n\in \mathbb N}` exists.
		4. *Left semi-infinite.* The backward orbit exists and does not repeat itself; but no part of the forward orbit exists.
		5. *Incomplete*. Only a finite number of images :math:`y\phi^{-n}, \dotsc, y\phi^{-1}, y, y\phi, \dotsc, y\phi^m` exist and all others do not.
	- a :py:func:`namedtuple <py3:collections.namedtuple>` consisting of
		1. type: one of the five types above
		2. power: if it exists, the characteristic power :math:`m` of this orbit; else ``None``.
		3. multiplier: if it exists, the characteristic multiplier :math:`\Gamme` of this orbit; else ``None``.
	
	 See definition 4.15 for details on characteristics.
	"""
	complete_infinite = 1
	complete_finite = 2
	right_semi_infinite = 3
	left_semi_infinite = 4
	incomplete = 5
	
	#Like vanilla namedtuples, do not allow Coord instances to have any instance attributes.
	__slots__ = ()
	
	#TODO docstrings
	@classmethod
	def periodic(cls, power):
		"""Shortcut for representing complete finite (i.e. periodic) orbits. Arguments are the characteristic of the orbit."""
		return cls(cls.complete_finite, power, None)
	
	@classmethod
	def semi_infinite(cls, direction, power, multiplier):
		"""Shortcut for representing periodic orbits. The argument *direction* should be one of the strings ``'left'`` or ``'right'``, indicating the direction in which this orbit is semi-infinite. The other two aArguments are the characteristic of the orbit.
		
		:raises ValueError: if direction is neither ``'left'`` nor ``'right'``.
		"""
		if direction == 'right':
			type = cls.right_semi_infinite
		elif direction == 'left':
			type = cls.left_semi_infinite
		else:
			raise ValueError("Must specify direction as either 'left' or 'right'.")
		return cls(type, power, multiplier)
	
	@classmethod
	def bi_infinite(cls):
		"""Shortcut for representing orbits which are infinite in both directions."""
		return cls(OrbitType.complete_infinite, None, None)

#todo renameme
def dump_orbit_types(aut, basis=None, words=None):
	r"""Prints the classification of the orbits under *aut* of each word in *words* with respect to *basis*. If *basis* is omitted, it is taken to be the minimal expansion given by :meth:`~thompson.automorphism._minimal_expansion`. If *words* is omited, it is taken to be the same as *basis*. See the docstring for :meth:`~thompson.automorphism._orbit_type`.
	
		>>> dump_orbit_types(arity_three_order_inf)
		x1 a1: OrbitType.left_semi_infinite
		x1 a2: OrbitType.complete_infinite
		x1 a3 a1: OrbitType.complete_infinite
		x1 a3 a2: OrbitType.complete_infinite
		x1 a3 a3: OrbitType.right_semi_infinite
		with respect to the basis [x1 a1, x1 a2, x1 a3 a1, x1 a3 a2, x1 a3 a3]
	"""
	if basis is None:
		basis = aut._minimal_expansion()
	if words is None:
		words = basis
	for w in words:
		print(w, ':', sep='', end=' ')
		print(aut._orbit_type(w, basis))
	print('with respect to the basis', basis)