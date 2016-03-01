"""
.. testsetup::
	
	from thompson.number_theory import gcd
	from thompson.word          import Signature, from_string
	from thompson.orbits        import print_component_types
	from thompson.automorphism  import *
	from thompson.examples      import *
"""
import os
import sys

from collections  import defaultdict
from copy         import copy
from fractions    import Fraction
from itertools    import chain, product

from .number_theory import lcm
from .word          import Signature, Word, free_monoid_on_alphas, format
from .generators    import Generators
from .homomorphism  import Homomorphism
from .orbits        import ComponentType, Characteristic, SolutionSet
from .pconj         import PowerCollection, search_pattern, mirrored
from .utilities     import handle_domain, intersection_from_domain
from .drawing.tpd   import generate_tikz_code

___all__ = ["Automorphism"]

class Automorphism(Homomorphism):
	r"""Represents an automorphism as a bijection between :meth:`bases <thompson.generators.Generators.is_basis>`.
	
	Generic attributes:
	
	:ivar signature: The :class:`~thompson.word.Signature` shared by the generating sets domain and range.
	:ivar quasinormal_basis: See :meth:`compute_quasinormal_basis`.
	:ivar pond_banks: A list of tuples :math:`(\ell, k, r)` such that :math:`(\ell, r)` are banks of a pond with :math:`\ell\phi^k = r`. For instance:
	
	.. doctest::
		
		>>> olga_f = load_example('olga_f')
		>>> olga_g = load_example('olga_g')
		>>> print(len(olga_f.pond_banks), sep=', ') #No ponds
		0
		>>> print(len(olga_g.pond_banks)) #One pond
		1
		>>> print(*olga_g.pond_banks[0], sep=',  ') 
		x1 a2 a1 a2 a1 a1,  2,  x1 a1 a2 a2
	
	Periodic attributes:
	
	:ivar periodic_orbits: a mapping :math:`d \mapsto L_d`, where :math:`L_d is the list of size :math:`d` orbits in the quasinormal basis.
	:ivar multiplicity: a mapping :math:`d \mapsto m_\phi(d, X_\phi)`, which is the number of periodic orbits of size :math:`d` in the :meth:`quasinormal basis <thompson.automorphism.Automorphism.compute_quasinormal_basis>` for :math:`\phi`. See Definition :paperref:`cycletypes` in the paper.
	:ivar cycle_type: the set :math:`\{d \in \mathbb{N} : \text{$\exists$ an orbit of length $d$.}\}`
	:ivar order: The :func:`~thompson.number_theory.lcm` of the automorphism's cycle type. This is the group-theoretic order of the :mod:`periodic factor <thompson.periodic>` of :math:`\phi`. If the cycle type is empty, the order is :math:`\infty`.
	
	..doctest::
		>>> def display_orbits(orbits_by_size):
		... 	for key in sorted(orbits_by_size):
		... 		print('Orbits of length', key)
		... 		for list in orbits_by_size[key]:
		... 			print('...', *list, sep=' -> ', end=' -> ...\n')
		>>> display_orbits(load_example('example_5_9').periodic_orbits)
		Orbits of length 2
		... -> x1 a2 a1 a1 -> x1 a2 a1 a2 -> ...
		... -> x1 a2 a2 a1 -> x1 a2 a2 a2 -> ...
		Orbits of length 3
		... -> x1 a1 a1 a1 -> x1 a1 a1 a2 -> x1 a1 a2 -> ...
	
	.. note::
		:mod:`mixed automorhpisms <thompson.mixed>` will have a **finite** order, despite being infinite-order group elements.
	
	Infinite attributes:
	
	:ivar characteristics: the set of characteristics :math:`(m, \Gamma)` of this automorphism.
	"""
	#Initialisation
	def __init__(self, domain, range, reduce=True):
		r"""Creates an automorphism mapping the given *domain* basis to the *range* basis in the given order.
		
		.. math:: \text{domain}_{\,i} \mapsto \text{range}_{\,i}
		
		The :meth:`quasi-normal basis <compute_quasinormal_basis>` :math:`X` and the various attributes are all calculated at creation time.
		
		:raises TypeError: if the bases have different arities or alphabet sizes.
		:raises TypeError: if either basis :meth:`isn't actually a basis <thompson.generators.Generators.is_basis>`.
		
		.. seealso:: :meth:`The superclass method <thompson.homomorphism.Homomorphism.__init__>`.
		"""
		#Check to see that the domain and range algebras are the same
		if domain.signature != range.signature:
			raise TypeError('Domain {} and range {} have different signatures.'.format(
			  domain.signature, range.signature))
		
		#Check to see that range is a basis for the range algebra
		#Have to expand non-simple words first
		Homomorphism._expand(domain, range)
		i, j = range.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		
		self._inv = {}
		super().__init__(domain, range, reduce)
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{}. Missing elements are {}.".format(
			  range.signature, [format(x) for x in missing]))
		
		#Everthing seems to be okay here.
		self.signature = domain.signature
		
		#Setup the inverse map
		for root in Generators.standard_basis(domain.signature):
			self._image_simple_above_domain(root, self.range.signature, self.domain.signature, self._inv)
		
		self.compute_quasinormal_basis()
		self.pond_banks = self._find_ponds()
		#This is extremely naughty.
		#Cast this class as a pure/mixed/infinite aut as appropriate.
		if len(self.multiplicity) == 0:
			from .infinite import InfiniteAut
			self.__class__ = InfiniteAut
		elif len(self.characteristics) == 0:
			from .periodic import PeriodicAut
			self.__class__ = PeriodicAut
		else:
			from .mixed import MixedAut
			self.__class__ = MixedAut
	
	#Computing images
	def _set_image(self, d, r, sig_in, sig_out, cache):
		if cache is not self._map and cache is not self._inv:
			raise ValueError("Incorrect cache provided.")
		super()._set_image(d, r, sig_in, sig_out, cache)
		cache = self._inv if cache is self._map else self._map
		cache[r] = d
	
	def image(self, key, inverse=False):
		"""If *inverse* is True, the inverse of the current automorphism is used to map *key* instead. Otherwise this method delegates to :meth:`Homomorphism.image <thompson.homomorphism.Homomorphism.image>`.
		
		Examples of finding inverse images:
		
			>>> phi = load_example('example_5_15')
			>>> print(phi.image('x1 a2 a2', inverse=True))
			x1 a2 a2 a1 a1
			>>> print(phi.image('x1 a1 a1 a2 a2 a1', inverse=True))
			x1 a2 a1 a2 a1
			>>> print(phi.image('x a2', inverse=True))
			x1 a2 a2 a2 x1 a2 a2 a1 a1 L
			>>> print(phi.image('x a2 a2 x a1 a2 L', inverse=True))
			x1 a2 a2 a1
		"""
		if inverse:
			return super().image(key, self.range.signature, self.domain.signature, self._inv)
		else:
			return super().image(key)
	
	def image_of_set(self, set, inverse=False):
		"""If *inverse* is True, the inverse of the current automorphism is used to map *set* instead. Otherwise this method delegates to :meth:`Homomorphism.image_of_set <thompson.homomorphism.Homomorphism.image_of_set>`.
		
		Examples of finding inverse images:
		
			>>> basis = Generators.standard_basis((2,1))
			>>> basis.expand_to_size(8);
			>>> print(load_example('example_5_3').image_of_set(basis, inverse=True))
			[x1 a1 a1 a1 a1, x1 a1 a1 a1 a2 x1 a1 a1 a2 L, x1 a1 a2 a2, x1 a1 a2 a1, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2 a1, x1 a2 a2 a2 a2]
		"""
		if inverse:
			return super().image_of_set(set, self.range.signature, self.domain.signature, self._inv)
		else:
			return super().image_of_set(set)
	
	def repeated_image(self, key, power):
		r"""If :math:`\psi` is the current automorphism, returns :math:`\text{key}\psi^\text{ power}`.
		
		:rtype: a :class:`~thompson.word.Word` instance.
		
		.. doctest::
			
			>>> phi = load_example('example_5_15')
			>>> print(phi.repeated_image('x1 a1', 10))
			x1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1
			>>> print(phi.repeated_image('x1 a1 a1 a1 a1 a1 a1 a1', -3))
			x1 a1
			>>> phi = load_example('arity_four')
			>>> print(phi.repeated_image('x1 a4 a4 a2', 4))
			x1 a3 a3 a2
			>>> print(phi.repeated_image('x1 a3 a3 a2', -4))
			x1 a4 a4 a2
		
		.. todo::
			This could be made more efficient for large values of *power* by using knowledge of the component containing *key*.
		"""
		inverse = power < 0
		if power == 0:
			if not isinstance(key, Word):
				key = Word(key, self.signature)
			return key
		power = abs(power)
		image = key
		for _ in range(power):
			image = self.image(image, inverse)
		return image
	
	#Group operations
	def __pow__(self, power):
		"""
		>>> psi = random_automorphism()
		>>> ~psi ** 2 == ~psi * ~psi == psi ** -2
		True
		"""
		if not isinstance(power, int):
			raise TypeError("Power must be an integer (instead of {}).".format(power))
		
		if power == 0:
			return Automorphism.identity(self.signature)
		if power == 1:
			return self
		if power == -1:
			return ~self
		
		domain = copy(self.domain)
		range = Generators(self.signature)
		for word in domain:
			range.append(self.repeated_image(word, abs(power)))
		if power > 0:
			pow = Automorphism(domain, range)
			pow.domain_relabeller = self.domain_relabeller
			pow.range_relabeller = self.range_relabeller
		else:
			pow = Automorphism(range, domain)
			pow.domain_relabeller = self.range_relabeller
			pow.range_relabeller  = self.domain_relabeller
		return pow
	
	def __invert__(self):
		"""We overload python's unary negation operator ``~`` as shorthand for inversion. (In Python, ``~`` is normally used for bitwise negation.) We can also call a method explicitily: ``phi.inverse()`` is exactly the same as ``~phi``.
		
			>>> phi = random_automorphism()
			>>> phi * ~phi == ~phi * phi
			True
			>>> (phi * ~phi).is_identity()
			True
			>>> (~phi).quasinormal_basis == phi.quasinormal_basis
			True
		"""
		inv = copy(self)
		inv.domain, inv.range = Generators.sort_mapping_pair(self.range, self.domain)
		inv._map, inv._inv = self._inv, self._map
		
		inv.pond_banks = {(r, -k, ell) for ell, k, r in self.pond_banks}
		inv.characteristics = {(-m, Gamma) for m, Gamma in self.characteristics}
		
		inv.domain_relabeller = self.range_relabeller
		inv.range_relabeller  = self.domain_relabeller
		return inv
	inverse = __invert__
	
	#Quasinormal basis
	#Some ideas for making this process faster:
	#- Cache the component types (wrt QNB) for the QNB elements
	#- Put some part of the core of type B orbits into confirmed? See /theory/Characteristic components and the QNB.tex
	#- Cache the powers (-ell, m) which determine the core part of a component
	def compute_quasinormal_basis(self):
		r"""We say that :math:`\phi` is *in semi-normal form* with respect to the basis :math:`X` if no element of :math:`X` lies in an incomplete :math:`X`-component of a :math:`\phi` orbit. See the :mod:`~thompson.orbits` module for more details.
		
		There is a minimal such basis, :math:`X_\phi` say, and we say that :math:`\phi` is *in quasi-normal form* with respect to :math:`X_\phi`. This method determines and returns the basis :math:`X_\phi` where :math:`\phi` denotes the current automorphism. The result is cached so that further calls to this method perform no additional computation. 
			
			>>> for name in ['example_4_5', 'alphabet_size_two', 'example_5_12_phi', 'example_6_2', 'example_6_8_phi']:
			... 	print(load_example(name).quasinormal_basis)
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			[x1 a1, x1 a2, x1 a3, x2]
			[x1 a1, x1 a2]
			[x1 a1 a1, x1 a1 a2, x1 a2]
			[x1 a1, x1 a2]
		
		:rtype: a :class:`~thompson.generators.Generators` instance.
		
		.. note:: This method is called automatically at creation time and is **not** needed to be called by the user.
		
		.. seealso:: Quasi-normal forms are introduced in Section :paperref:`sec:qnf` of the paper. In particular, this method implements Lemma :paperref:`lem:qnf`. Higman first described the idea of quasi-normal forms in Section 9 of [Hig74]_.
		"""
		#todo:: Make the QNB read-only somehow, so that it cannot be expanded once computed.
		self.quasinormal_basis = None
		self.pond_banks = None
		self.periodic_orbits = defaultdict(list)
		self.multiplicity = {}
		
		#1. Expand the starting basis until each no element's belongs to a finite X-component.
		basis = self.seminormal_form_start_point()
		basis.cache = set(basis)
		confirmed = set()
		i = 0
		checks_needed = len(basis)
		
		while checks_needed > 0:
			if basis[i] in confirmed:
				i = (i + 1) % len(basis)
				checks_needed -= 1
				continue
			
			ctype, images, _ = self.orbit_type(basis[i], basis)
			if ctype.is_incomplete():
				basis._expand_with_cache(i)
				checks_needed = len(basis)
				continue
			
			elif ctype.is_type_A():
				confirmed.update(images.values())
				#expand basis until each of the images is below basis
				for img in images.values():
					#TODO. I think this operation takes O(len(basis)).
					index, tail = basis.test_above(img, return_index=True)
					for j in range(len(tail)):
						basis._expand_with_cache(index)
						index += -tail[j] - 1
						checks_needed += self.signature.arity - 1
					assert img in basis
				
				#record this new periodic orbit
				period = ctype.characteristic[0]
				orbit = tuple( images[i] for i in range(period) )
				self.periodic_orbits[period].append(orbit)
			
			elif ctype.is_type_B():
				confirmed.add(basis[i])
			i = (i + 1) % len(basis)
			checks_needed -= 1
		
		#Tidy up data
		self.quasinormal_basis = basis
		self.multiplicity = {
			period: len(orbits)
			for period, orbits in self.periodic_orbits.items()
		}
		self.cycle_type = frozenset(self.multiplicity)
		self.order = float('inf') if len(self.cycle_type) == 0 else lcm(self.cycle_type)
		
		self.characteristics = set()
		for endpt in chain(*self.semi_infinite_end_points()):
			ctype, _, _ = self.orbit_type(endpt, basis)
			if ctype.is_type_B():
				self.characteristics.add(ctype.characteristic)
	
	def _find_ponds(self):
		terminal, initial = self.semi_infinite_end_points(exclude_characteristics=True)
		terminal_data = set()
		for term in terminal:
			tail = self._descend_to_complete_infinite(term)
			terminal_data.add((term, tail))
		initial = set(initial)
		
		ponds = []
		for (term, tail) in terminal_data:
			for init in initial:
				power = self._are_banks(term, init, tail)
				if power is not None:
					ponds.append((term, power, init))
					initial.remove(init)
					break
		return ponds
	
	#TODO: maybe this ought to be part a method of Generators or Signature
	def _intersection_of_trees(self, domain, range):
		r"""Given the leaves A and B of two trees (i.e. bases), computes the leaves (basis) X of the tree intersection A & B.
		"""
		basis = Generators.standard_basis(self.signature)
		i = 0
		while i < len(basis):
			b = basis[i]
			if b in domain or b in range:
				i += 1
			else:
				basis.expand(i)
		return basis
	
	def seminormal_form_start_point(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \cup Z`. This basis that this method produces is the smallest possible which *might* be semi-normal.
		
			>>> for name in ['example_4_5', 'example_4_11', 'example_4_12', 'example_5_15', 'cyclic_order_six']:
			... 	print(load_example(name).seminormal_form_start_point())
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			[x1 a1, x1 a2]
			[x1 a1, x1 a2]
			[x1 a1, x1 a2 a1, x1 a2 a2]
			[x1 a1 a1, x1 a1 a2 a1, x1 a1 a2 a2, x1 a2]
		
		.. seealso::  Remark :paperref:`rem:snf_start_point` of the paper.
		"""
		return self._intersection_of_trees(self.domain, self.range)
	
	def orbit_type(self, y, basis=None):
		"""Returns the orbit type of *y* with respect to the given *basis*. If *basis* is omitted, the :meth:`quasi-normal basis <compute_quasinormal_basis>` is used by default. Also returns a dictionary of computed images, the list (:paperref:`eq:uorb`) from the paper.
		
			>>> phi = load_example('example_4_5')
			>>> print_component_types(phi, phi.domain)
			x1 a1 a1 a1: Left semi-infinite component with characteristic (-1, a1)
			x1 a1 a1 a2: Bi-infinite component
			x1 a1 a2: Right semi-infinite component with characteristic (1, a2)
			x1 a2 a1: Periodic component of order 2
			x1 a2 a2: Periodic component of order 2
			with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> print_component_types(phi, basis=phi.domain, words=['x', 'x a1', 'x a2'])
			x1: Incomplete component
			x1 a1: Incomplete component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			
		.. doctest::
			:hide:
			
			>>> print_component_types(load_example('example_4_11'))
			x1 a1: Left semi-infinite component with characteristic (-1, a1)
			x1 a2: Right semi-infinite component with characteristic (1, a2)
			with respect to the basis [x1 a1, x1 a2]
		
		.. doctest::
			:hide:
			
			>>> phi = load_example('example_4_12')
			>>> basis = phi.seminormal_form_start_point()
			>>> print_component_types(phi, basis)
			x1 a1: Incomplete component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1, x1 a2]
			>>> print(basis.expand(0))
			[x1 a1 a1, x1 a1 a2, x1 a2]
			>>> print_component_types(phi, basis)
			x1 a1 a1: Bi-infinite component
			x1 a1 a2: Bi-infinite component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2]
			>>> print(basis.expand(2))
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> print_component_types(phi, basis)
			x1 a1 a1: Periodic component of order 4
			x1 a1 a2: Periodic component of order 4
			x1 a2 a1: Periodic component of order 4
			x1 a2 a2: Periodic component of order 4
			with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> #Example 4.25
			>>> print_component_types(load_example('example_5_15'))
			x1 a1: Right semi-infinite component with characteristic (1, a1 a1)
			x1 a2 a1: Bi-infinite component
			x1 a2 a2: Left semi-infinite component with characteristic (-1, a1 a1)
			with respect to the basis [x1 a1, x1 a2 a1, x1 a2 a2]
		
		Components can be calculated with respect to any *basis*, not just the :meth:`quasi-normal basis <compute_quasinormal_basis>`.
		
			>>> print(olga_f.quasinormal_basis)
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2]
			>>> basis = olga_f.quasinormal_basis.copy().expand(-1); print(basis)
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2 a1, x1 a2 a2 a2 a2]
			>>> from pprint import pprint
			>>> print(olga_f.orbit_type('x a2 a2 a2')[0])
			Bi-infinite component
			>>> print(olga_f.orbit_type('x a2 a2 a2', basis)[0])
			Incomplete component
		
		.. doctest::
			:hide:
			
			>>> ctype, images, _ = olga_f.orbit_type('x a1 a2 a1')
			>>> print(ctype); pprint(images)
			Left semi-infinite component with characteristic (-1, a1)
			{-1: Word('x1 a1 a2 a1 a1', (2, 1)),
			 0: Word('x1 a1 a2 a1', (2, 1)),
			 1: Word('x1 a1 a2', (2, 1))}
			>>> ctype, images, _ = olga_f.orbit_type('x a1 a2 a1', basis)
			>>> #Component of word not changed
			>>> print(ctype); pprint(images)
			Left semi-infinite component with characteristic (-1, a1)
			{-1: Word('x1 a1 a2 a1 a1', (2, 1)),
			 0: Word('x1 a1 a2 a1', (2, 1)),
			 1: Word('x1 a1 a2', (2, 1))}
		
		.. seealso:: Lemmas :paperref:`ABC`, :paperref:`lem:qnf` of the paper.
		
		.. todo:: This should be renamed to component_type.
		"""
		images = {}
		type_b_data = None
		
		if isinstance(y, str):
			y = Word(y, (2, 1))
		
		if basis is None and self.quasinormal_basis is not None:
			basis = self.quasinormal_basis
		rinf, rpow1, rpow2, rimages = self.test_semi_infinite(y, basis, backward=False)
		for i, image in enumerate(rimages):
			images[i] = image
		
		#Eliminate the periodic case first.
		if rinf and rimages[rpow2] == y:
			assert rpow1 == 0
			ctype = ComponentType.periodic(rpow2)
			del images[rpow2]
			return ctype, images, type_b_data
		
		#Else, we don't have anything periodic.
		linf, lpow1, lpow2, limages = self.test_semi_infinite(y, basis, backward=True)
		for i, image in enumerate(limages):
			images[-i] = image
		
		if not (linf or rinf):
			ctype = ComponentType.incomplete()
			
		elif linf and not rinf:
			tail = y.test_above(limages[-1])
			if tail is None:
				characteristic = None
				type_b, tail = basis.test_above(limages[lpow1])
				type_b_data = type_b_triple(-lpow1, type_b, tail)
			else:
				assert len(tail) > 0
				assert lpow1 == 0
				characteristic = Characteristic(-lpow2, tail)
			ctype = ComponentType.semi_infinite(characteristic, backward=True)
		
		elif rinf and not linf:
			tail = y.test_above(rimages[-1])
			if tail is None:
				characteristic = None
				type_b, tail = basis.test_above(rimages[rpow1])
				type_b_data = type_b_triple(rpow1, type_b, tail)
			else:
				assert len(tail) > 0
				assert rpow1 == 0
				characteristic = Characteristic(rpow2, tail)
			ctype = ComponentType.semi_infinite(characteristic, backward=False)
		
		elif linf and rinf:
			ctype = ComponentType.complete_infinite()
			type_b, tail = basis.test_above(rimages[rpow1])
			type_b_data = type_b_triple(rpow1, type_b, tail)
		
		#Test_semi_infinite returns one more word than is needed for the 'core'.
		#When appropriate, we remove this extra word.
		#TODO this breaks the conjugacy test and needs to be implemented properly.
		# for index in (rpow2, -lpow2):
			# if abs(index) > 0:
				# try:
					# del images[index]
				# except KeyError:
					# pass
		return ctype, images, type_b_data
	
	def test_semi_infinite(self, y, basis, backward=False):
		r"""Computes the orbit type of *y* under the current automorphism :math:`\psi` with respect to *basis* in the given direction. Let :math:`y\psi^m` be the most recently computed image. The process stops when either:
		
			1. :math:`y\psi^m` is not below the *basis*, for some :math:`m\geq 0`.
			
				- infinite: ``False``
				- start: ``0``
				- end: ``m-1``
				- images: :math:`y, y\psi, \dotsc, y\psi^{m-1}`.
			
			2. For some :math:`0 \le l < m`, :math:`y\psi^l` is above (or equal to) :math:`y\psi^m`.
			
				- infinite: ``True``
				- start: ``l``
				- end: ``m``
				- images: :math:`y, y\psi, \dotsc, y\psi^{m}`.
				
			.. note::
				
				The word :math:`y\psi^m` is not strictly in the core part of the orbit of Lemma 4.24. We return this as part of *images* so that we can compute the characteristic multiplier in :meth:`orbit_type`.
		
		:returns: the tuple *(infinite, start, end, images)*.
		
		.. seealso:: Lemma :paperref:`lem:qnf` of the paper.
		"""
		image = y
		images = [y]
		m = 0
		heads = set()
		
		result = basis.test_above(y)
		if result is None:
			return False, 0, m-1, images 
		
		prefix, _ = result
		while True:
			heads.add(prefix)
			m += 1
			image = self.image(image, inverse=backward) #Compute y\phi^{\pm m}
			result = basis.test_above(image)
			#1. Did we fall out of X<A>?
			if result is None:
				return False, 0, m-1, images
			
			#2. Otherwise, is there evidence to conclude that this is this an infinite orbit?
			prefix, tail = result
			if prefix in heads:
				for ell, previous in enumerate(images):
					if len(image) >= len(previous) and prefix == previous[:len(prefix)]:
						images.append(image)
						return True, ell, m, images
			images.append(image)
	
	def semi_infinite_end_points(self, exclude_characteristics=False):
		r"""Returns the list of terminal :class:`Words <thompson.word.Word>` in left semi-infinite components and the list of initial words in right semi-infinite components. This is all computed with respect the current automorphism's :meth:`quasinormal basis <compute_quasinormal_basis>`. These are the sets :math:`X\langle A\rangle \setminus Y\langle A\rangle` and :math:`X\langle A\rangle \setminus Z\langle A\rangle`.
		
		:param exclude_characteristics: if True, only the endpoints of non-characteristic semi-infinite orbits will be returned.
		
		.. doctest::
			
			>>> for id in ['4_5', '4_11', '4_12', '5_15', '6_2']:
			... 	aut = load_example('example_' + id)
			... 	print(*aut.semi_infinite_end_points())
			[x1 a1 a1] [x1 a1 a2]
			[x1 a1] [x1 a2]
			[] []
			[x1 a2 a2, x1 a2 a2 a1] [x1 a1, x1 a1 a1]
			[x1 a1 a1] [x1 a2]
			>>> phi = load_example('nathan_pond_example')
			>>> print(*phi.semi_infinite_end_points())
			[x1 a1 a1, x1 a1 a1 a1, x1 a1 a1 a2] [x1 a1 a2, x1 a1 a2 a1, x1 a1 a2 a2]
			>>> print(*phi.semi_infinite_end_points(exclude_characteristics=True))
			[x1 a1 a1 a2] [x1 a1 a2 a2]
		
		:rtype: A pair of :class:`Generators <thompson.generators.Generators>`.
		
		.. seealso:: The discussion before Lemma :paperref:`9.1H`.
		"""
		basis = self.quasinormal_basis                    #X
		min_expansion = basis.minimal_expansion_for(self) #Y
		img_expansion = self.image_of_set(min_expansion)  #Z
		terminal = basis.descendants_above(min_expansion) #X<A> \ Y<A>
		initial  = basis.descendants_above(img_expansion) #X<A> \ Z<A>
		
		if exclude_characteristics:
			terminal = Generators(self.signature, (t for t in terminal if not self.orbit_type(t, basis)[0].is_type_B()))
			initial =  Generators(self.signature, (i for i in initial  if not self.orbit_type(i, basis)[0].is_type_B()))
		return terminal, initial
	
	
	def _descend_to_complete_infinite(self, endpt):
		"""A really crude use of the lemma described in AJD's email.
		
		**Lemma.** Let :math:`b` be the bank of a pond. There is a :math:`\Gamma \in A^*` such that :math:`b\Gamma` belongs to a doubly-infinite :math:`\phi`-orbit that meets :math:`X`.
		"""
		basis = self.quasinormal_basis
		for gamma in free_monoid_on_alphas(self.signature.arity):
			w = endpt.extend(gamma)
			ctype, _, _ = self.orbit_type(w, basis)
			if ctype.type == ComponentType._complete_infinite:
				return gamma
	
	def _are_banks(self, term, init, tail):
		#todo docstring
		u = term.extend(tail)
		v = init.extend(tail)
		solns = self.share_orbit(u, v)
		assert not solns.is_sequence(), solns
		if solns.is_empty():
			return None
		power = solns.base
		if self.repeated_image(term, power) == init:
			return power
		else:
			return None
	
	def dump_QNB(self):
		"""A convenience method for printing the quasinormal basis :math:`X`. The :math:`X`-components of the elements :math:`x \in X` are displayed.
		
			>>> load_example('example_5_3').dump_QNB()
			x1 a1 a1 a1 Left semi-infinite component with characteristic (-1, a1)
			x1 a1 a1 a2 Right semi-infinite component with characteristic (1, a2)
			x1 a1 a2 a1 Periodic component of order 2
			x1 a1 a2 a2 Periodic component of order 2
			x1 a2 a1 Right semi-infinite component with characteristic (1, a1)
			x1 a2 a2 Left semi-infinite component with characteristic (-1, a2)
		"""
		for gen in self.quasinormal_basis:
			ctype, _, _ = self.orbit_type(gen, self.quasinormal_basis)
			print(gen, ctype)
	
	def print_characteristics(self):
		"""For convenience, a method that prints out all of the characteristics of type B components wrt the quasinormal_basis.
		
		.. doctest::
			:options: +NORMALIZE_WHITESPACE
			
			>>> for id in ['4_1', '5_15', '6_2', '6_8_phi']:
			... 	name = 'example_' + id
			... 	print(name)
			... 	load_example(name).print_characteristics()
			example_4_1
			(-1, a1)
			(1, a2)
			example_5_15
			(-1, a1 a1)
			(1, a1 a1)
			example_6_2
			(-2, a1)
			(1, a2)
			example_6_8_phi
			(-1, a1 a1 a1)
			(1, a2 a2 a2)
			>>> (load_example('example_6_2')**2).print_characteristics()
			(-1, a1)
			(1, a2 a2)
			>>> #Lemma 5.16
			>>> psi, phi = random_conjugate_pair()
			>>> psi.characteristics == phi.characteristics
			True
		
		.. doctest::
			:hide:
			
			>>> #Lemma 6.1
			>>> from random import randint
			>>> psi = random_infinite_automorphism()
			>>> original_chars = psi.characteristics
			>>> power = psi
			>>> a = randint(2, 6)
			>>> for _ in range(a - 1): 
			... 	power *= psi
			>>> chars = set()
			>>> for mult, char in original_chars:
			... 	d = gcd(mult, a)
			... 	q = abs(a // d)
			... 	chars.add((mult//d, char*q))
			>>> chars == power.characteristics
			True
		
		.. seealso:: Defintion :paperref:`setMultipliers`.
		"""
		for power, mult in sorted(self.characteristics):
			print('({}, {})'.format(power, format(mult)))
	
	#Orbit sharing test
	def share_orbit(self, u, v):
		r"""Determines if :math:`u` and :math:`v` are in the same orbit of the current automorphism :math:`\psi`. Specifically, does there exist an integer :math:`m` such that :math:`u\psi^m = v`?
		
		.. doctest::
			:options: -ELLIPSIS
			
			>>> phi = load_example('alphabet_size_two')
			>>> u  = Word('x1 a2 a3 a1 a2', (3, 2))
			>>> v1 = Word('x1 a1 a2 a2 a3 a1', (3, 2))
			>>> v2 = Word('x2 a3 a2', (3, 2))
			>>> print(phi.share_orbit(u, v1))
			{}
			>>> print(phi.share_orbit(u, v2))
			{-2}
			>>> print(phi.share_orbit(u, u))
			{0}
		
		.. doctest::
			:hide:
			:options: -ELLIPSIS
			
			>>> phi = load_example('example_5_15')
			>>> u  = Word('x a2 a2 a1 a1 a2', (2, 1))
			>>> v1 = Word('x a1 a2', (2, 1))
			>>> v2 = Word('x a1 a1 a2', (2, 1))
			>>> print(phi.share_orbit(u, v1))
			{}
			>>> print(phi.share_orbit(u, v2))
			{3}
			>>> u  = Word('x a2 a2 x a1 a2 L x a2 a1 L x a1 a1 a1 a2 L', (2, 1))
			>>> vs = [
			... 	Word('x a2 a2 a1 a1 a1 a1 a1', (2, 1)),
			... 	Word('x a2 a2 a2', (2, 1)),
			... 	Word('x a2 a2 x a1 a2 L', (2, 1)),
			... 	Word('x a1 a1 x a1 a2 x a2 a2 L L', (2, 1))]
			... 
			>>> for v in vs: print(phi.share_orbit(u, v))
			{-4}
			{}
			{-1}
			{}
			>>> phi = load_example('arity_four')
			>>> u1 = Word('x a2 a3 a1', (4, 1))
			>>> v1 = Word('x a3 a3 a3 a3 a3 a3 a3 a3 a3 a2 a3 a1', (4, 1))
			>>> v2 = Word('x a1 a2 a3 a4', (4, 1))
			>>> u2 = Word('x a3 a4 a1', (4, 1))
			>>> v3 = Word('x a4 a1 a1', (4, 1))
			>>> v4 = Word('x a4 a3 a2 a1', (4, 1))
			>>> print(phi.share_orbit(u1, v1))
			{9}
			>>> print(phi.share_orbit(u1, v2))
			{}
			>>> print(phi.share_orbit(u2, v3))
			{-1}
			>>> print(phi.share_orbit(u2, v4))
			{}
			>>> u = Word("x a1 a2 a1 a2 a1", (2, 1))
			>>> v = Word("x a2 a2 a2 a1", (2, 1))
			>>> phi = load_example('cyclic_order_six')
			>>> print(phi.share_orbit(u, v))
			{..., -1, 2, 5, 8, 11, 14, ...}
			>>> u = Word("x a1 a1 x a1 a2 a1 x a1 a2 a2 a1 L L", (2, 1))
			>>> v1 = Word("x a1 a1 x a2 a2 x a2 a1 L L", (2, 1))
			>>> v2 = Word("x a1 a1 x a1 a2 a1 x a1 a1 L L", (2, 1))
			>>> print(phi.share_orbit(u, v1))
			{..., -1, 5, 11, 17, 23, 29, ...}
			>>> print(phi.share_orbit(u, v2))
			{}
			>>> #Two sides of a pond
			>>> u = Word('x a1 a1 a1 a1 a1 a1 a1 a2', (2, 1))
			>>> v = Word('x a2 a2 a1 a1 a2', (2, 1))
			>>> print(load_example('first_pond_example_phi').share_orbit(u, v))
			{4}
		
		:returns: The (possibly empty) :class:`~thompson.orbits.SolutionSet` of all integers :math:`m` for which :math:`u\psi^m = v`. Note that if :math:`u = v` this method returns :math:`\mathbb{Z}`. 
		
		.. seealso:: The implementation is due to lemma :paperref:`9.7H` of the paper.
		"""
		#TODO a script which randomly checks examples to verify.
		orig_u, orig_v = u, v
		basis = self.quasinormal_basis
		if not (basis.is_above(u) and basis.is_above(v)):
			depth = max(u.max_depth_to(basis), v.max_depth_to(basis))
			alphas = range(-1, -self.signature.arity - 1, -1)
			solution_set = SolutionSet.the_integers()
			# For all strings \Gamma of length *depth*:
			for tail in product(alphas, repeat=depth):
				u_desc = u.extend(tail)
				v_desc = v.extend(tail)
				solution_set &= self.share_orbit(u_desc, v_desc)
				if solution_set.is_empty():
					return solution_set
			return solution_set
			#A recursive alternative: intersect and return the results of tests on only the children of u, v.
		
		#Now we're dealing with simple words below the basis.
		u_head, u_tail = basis.test_above(u)
		v_head, v_tail = basis.test_above(v)
		u_head_type, _, u_head_data = self.orbit_type(u_head, basis)
		v_head_type, _, v_head_data = self.orbit_type(v_head, basis)
		
		#Is either element periodic?
		if u_head_type.is_type_A() or v_head_type.is_type_A():
			#If so, are they both peroidic with the same period?
			if u_head_type != v_head_type:
				return SolutionSet.empty_set()
			
			period = u_head_type.characteristic[0]
			image = u
			for i in range(1, period+1):
				image = self.image(image)
				if image == v:
					return SolutionSet(i, period)
			return SolutionSet.empty_set()
		
		#Neither u nor v are periodic.
		#Move along the orbit until we find a nicer u and v.
		u_shift, u_head, u_tail = self._type_b_descendant(u_head, u_tail, u_head_data)
		v_shift, v_head, v_tail = self._type_b_descendant(v_head, v_tail, v_head_data)
		
		u = u_head.extend(u_tail)
		v = v_head.extend(v_tail)
		type, images, _ = self.orbit_type(u, basis)
		
		if self.pond_banks is not None:
			#Add on the other side of a pond if neccessary
			self._check_for_ponds(images)
		for i, image in images.items():
			if image == v:
				assert self.repeated_image(orig_u, u_shift + i - v_shift) == orig_v, (format(orig_u), u_shift + i - v_shift)
				return SolutionSet.singleton(u_shift + i - v_shift)
		return SolutionSet.empty_set()
	
	def _type_b_descendant(self, head, tail, head_data):
		r"""Takes a pair :math:`(y, \Gamma)` below the quasi-normal basis :math:`X` and returns a triple :math:`(\widetilde{y}, \widetilde{\Gamma}, k) where
		* The orbit of :math:`\widetilde{y}` is of type B;
		* :math:`\widetilde{y}\widetilde{\Gamma}` is in the same orbit as :math:`y\Gamma`.
		* :math:`\widetilde{\Gamma}` does not start with the characteristic multiplier of :math:`\widetilde{y}`.
		
			>>> phi = load_example('example_5_15')
			>>> basis = phi.quasinormal_basis
			>>> print(basis)
			[x1 a1, x1 a2 a1, x1 a2 a2]
			>>> head = Word('x a2 a2', (2, 1))
			>>> _, _, type_b_data = phi.orbit_type(head, basis)
			>>> phi._type_b_descendant(head, from_string('a1 a1 a2'), type_b_data)
			(1, Word('x1 a2 a2', (2, 1)), (-2,))
		"""
		orig_head = head
		orig = head.extend(tail)
		if head_data is not None:
			shift_1, head, prefix = head_data['power'], head_data['target'], head_data['end_tail']
			assert self.repeated_image(orig_head, shift_1) == head.extend(prefix), (orig_head, shift_1, head, prefix)
			tail = prefix + tail
		else:
			shift_1 = 0
		
		head_type, _, _ = self.orbit_type(head, self.quasinormal_basis)
		assert head_type.is_type_B(), (head, head_type)
		shift_2, tail = self._remove_from_start(tail, head_type.characteristic)
		assert self.repeated_image(orig, shift_1 + shift_2) == head.extend(tail), (shift_1, shift_2)
		return (shift_1 + shift_2), head, tail
	
	@staticmethod
	def _remove_from_start(tail, characteristic):
		"""Makes a copy of *tail* and removes as many copies of *multiplier* from the start of the copy as possible. Returns a tuple *(k, tail copy)* where *k* is the number of steps forward into the orbit we made.
			
			>>> tail = from_string('a1 a2 a1 a2 a1 a2 a3')
			>>> multiplier = from_string('a1 a2')
			>>> Automorphism._remove_from_start(tail, (2, multiplier))
			(-6, (-3,))
			>>> multiplier = from_string('a1')
			>>> Automorphism._remove_from_start(tail, (2, multiplier))
			(-2, (-2, -1, -2, -1, -2, -3))
		"""
		power, multiplier = characteristic
		i = 0
		n = len(multiplier)
		while tail[i*n : i*n + n] == multiplier:
			i += 1
		return -power*i, tail[i*n:]
	
	def _check_for_ponds(self, images):
		#todo doctests
		M = max(images.keys())
		m = min(images.keys())
		basis = self.quasinormal_basis
		
		for (ell, k, r) in self.pond_banks:
			if images[M] == ell:
				core = self.test_semi_infinite(r, basis)[3]
				for i, img in enumerate(core):
					images[M + k + i] = img
				return
			
			elif images[m] == r:
				core = self.test_semi_infinite(ell, basis, backward=True)[3]
				for i, img in enumerate(core):
					images[m - k - i] = img
				return
	
	#Other functions useful for ordinary conjugacy
	def is_conjugate_to(self, other):
		"""A shortcut for ``self.test_conjugate_to(other) is not None``."""
		return self.test_conjugate_to(other) is not None
	
	#Power conjugacy
	def _test_power_conjugate_upto(self, other, sbound, obound, inverses=False, cheat=False):
		r"""In both the periodic and infinite cases, we establish bounds on the powers :math:`a, b` for conjugacy; the rest is brute force.  This method tests to see if :math:`\psi^a` is conjugate to :math:`\phi^b` within the supplied bounds. Should it find a conjugator :math:`\rho`, this method yields a triple :math:`(a, b, \rho)`.
		
		Let :math:`\hat a, \hat b` denote *sbound* and *obound* respectively. If *inverses* is False, then we search over the ranges :math:`1 \le a \le \hat a` and :math:`1 \le b \le \hat b`. If *inverses* is True, we search over the (four times larger) range `1 \le |a| \le \hat a` and :math:`1 \leq |b| \le \hat b`.
		"""
		if cheat:
			from thompson.examples.random import random_power_bounds
			sbound, obound = random_power_bounds
		if sbound < 0:
			raise ValueError('sbound parameter should be at least 0 (received {}).'.format(sbound))
		if obound < 0:
			raise ValueError('obound parameter should be at least 0 (received {}).'.format(obound))
		# print(self.__class__.__name__, other.__class__.__name__, sbound, obound)
		if sbound == 0 or obound == 0:
			raise StopIteration
		
		#WLOG let sbound <= obound
		if sbound > obound:
			self, other = other, self
			sbound, obound = obound, sbound
			# print('swapping bounds to', sbound, obound)
			swapped = True
		else:
			swapped = False
		
		s_powers = PowerCollection(self)
		o_powers = PowerCollection(other)
		
		iterator = search_pattern(sbound, obound)
		if inverses:
			iterator = mirrored(iterator)
		for a, b in iterator:
			# print('trying', a, b)
			rho = s_powers[a].test_conjugate_to(o_powers[b])
			if rho is not None:
				yield (a, b, rho) if not swapped else (b, a, ~rho)
				if inverses:
					yield (-a, -b, rho) if not swapped else (-b, -a, ~rho)
	
	def preserves_order(self):
		"""Returns True if this automorphism is an element of :math:`F_{n,r}`, Higman's analogue of Thompson's group :math:`F`. Otherwise returns False.
		
			>>> random_automorphism(group='F').preserves_order()
			True
			>>> phi = random_automorphism(group='T')
			>>> #phi preserves order iff it is in F.
			>>> (sorted(phi.range) == phi.range) == phi.preserves_order()
			True
			>>> load_example('nathan_pond_example').preserves_order()
			False
		"""
		indices = range(len(self.range) - 1)
		return all(self.range[i] < self.range[i+1] for i in indices)
	
	def cycles_order(self):
		"""Returns True if this automorphism is an element of :math:`T_{n,r}`, Higman's analogue of Thompson's group :math:`T`. Otherwise returns False.
		
		>>> random_automorphism(group='F').cycles_order()
		True
		>>> random_automorphism(group='T').cycles_order()
		True
		>>> load_example('nathan_pond_example').cycles_order() # in V \ T
		False
		"""
		indices = range(len(self.range) - 1)
		looped = False
		for i in indices:
			if not self.range[i] < self.range[i+1]:
				if not looped and all(abs(x) == 1 for x in self.range[i+1]):
					looped = True
				else:
					return False
		return True
	
	def _end_of_iac(self, root, leaves, backward=False):
		"""Given a root :math:`r` of :math:`A \setminus B` of :math:`B \setminus A`, this finds the IAC containing :math:`r` and returns the endpoint :math:`u_1` or :math:`f(u_n)` which is not :math:`r`."""
		assert root in leaves
		u = root
		while u in leaves:
			u = self.image(u, inverse=backward)
		return u
	
	def test_revealing(self, domain='wrt QNB'):
		r"""Determines if the given automorphism :math:`\phi` is revealed by the tree pair :math:`(D, \phi(D))`, where :math:`D` is the given *domain*.
		
		The *domain* may be implicitly specified by the string ``'minimal'`` or ``'wrt QNB'``. In the first case, *domain* is taken to be the minimal *domain* required to specify the automorphism. In the second case, *domain* is taken to be the minimal expansion of the quasinormal basis.
		
		:returns: `None` if the pair is revealing for :math:`\phi`. Otherwise, returns (as a :class:`~thompson.word.Word`) the root of a component of either :math:`D \setminus \phi(D)` or :math:`\phi(D) \setminus D` which does not contain an attractor/repeller.
		
		>>> load_example('olga_f').test_revealing() is None
		True
		>>> print(load_example('semi_inf_c').test_revealing())
		x1 a2 a1
		>>> print(load_example('non_revealing').test_revealing())
		x1 a1 a1
		>>> f = load_example('cyclic_order_six')
		>>> print(f.test_revealing('minimal'))
		x1 a2
		>>> f.test_revealing('wrt QNB') is None
		True
		
		.. caution:: This is an experimental feature based on [SD10]_.
		"""
		domain = handle_domain(domain, self)
		range, X = intersection_from_domain(domain, self)
		
		roots_d_minus_r = X.filter(lambda x: x not in domain)
		roots_r_minus_d = X.filter(lambda x: x not in range )
		
		for roots, leaves in ( (roots_d_minus_r, range), (roots_r_minus_d, domain) ):
			for root in roots:
				w = self._end_of_iac(root, leaves, backward=( leaves == range ))
				if not root.is_above(w):
					return root
		return None
	
	def is_revealing(self, domain='wrt QNB'):
		r"""Calls :meth:`test_revealing`, but only returns ``True`` or ``False``.
		
		:returns: True if the pair is revealing, otherwise False.
		
		>>> load_example('olga_f').is_revealing()
		True
		>>> load_example('semi_inf_c').is_revealing()
		False
		>>> load_example('non_revealing').is_revealing()
		False
		>>> f = load_example('cyclic_order_six')
		>>> f.is_revealing('minimal')
		False
		>>> f.is_revealing('wrt QNB')
		True
		
		.. caution:: This is an experimental feature based on [SD10]_.
		"""
		return self.test_revealing(domain) is None
	
	def area_to_identity(self, scaled=False):
		r"""Let :math:`\phi \in F_{n, 1}`, viewed as a bijection of the interval.
		What is the (unsigned) area between :math:`\phi` and the identity?
		
		:param bool scaled: If true, the area is normalised by ignoring subintervals where the function is equal to the identity.
		
		.. warning:: This is an experimental feature based on a suggestion by Henry Bradford.
		
		.. doctest::
			
			>>> for n in range(4):
			... 	print(standard_generator(n).area_to_identity())
			... 
			5/32
			5/128
			5/512
			5/2048
			>>> x0 = standard_generator(0)
			>>> x1 = standard_generator(1)
			>>> g = ~x0 * x1 * x0
			>>> print(g.area_to_identity())
			3/32
			>>> f = load_example('non_dyadic_fixed_point')
			>>> print(f.area_to_identity())
			43/768
		
			>>> for n in range(4):
			...		print(standard_generator(n).area_to_identity(scaled=True))
			5/32
			5/32
			5/32
			5/32
		"""
		if not self.preserves_order():
			raise ValueError('This automorphism is not in F')
		
		#1. Work out where the function is different from the identity.
		num_sections = len(self.domain)
		non_identity_sections = []
		indices = []
		for i in range(num_sections):
			x0, x1 = self.domain[i].as_interval()
			y0, y1 = self.range[i].as_interval()
			matches_identity = (x0 == y0 and x1 == y1)
			if not matches_identity:
				indices.append(i)
			
			if matches_identity or i == num_sections - 1:
				if indices:
					non_identity_sections.append(indices)
					indices = []
		
		if not non_identity_sections:
			return Fraction(0, 1)
		
		#2. Compute the width and area of each non-identity section
		widths = []
		areas = []
		for indices in non_identity_sections:
			start = self.domain[indices[0]].as_interval()[0]
			end   = self.domain[indices[-1]].as_interval()[1]
			widths.append(end-start)
			area = 0
			for i in indices:
				x0, x1 = self.domain[i].as_interval()
				y0, y1 = self.range[i].as_interval()
				#Does this interval cross the diagonal?
				upward_crossing   = y0 < x0 and y1 > x1
				downward_crossing = y0 > x0 and y1 < x1
				if upward_crossing or downward_crossing:
					#compute the crossing point
					gradient = (y1 - y0) / (x1 - x0)
					intercept = y0 - gradient * x0
					intersection = intercept / (1 - gradient)
					#add the area of the two triangles
					area += trapezium_area(x0, y0, intersection, intersection)
					area += trapezium_area(intersection, intersection, x1, y1)
				else:
					area += trapezium_area(x0, y0, x1, y1)
			areas.append(area)
		
		if not scaled:
			scalar = 1
		else:
			scalar = sum(widths) ** -2
		return sum(areas) * scalar
	
	def centralise_period(self, period, trees, rearranger):
		"""Constructs a new element :math:`\psi` commuting with the given element :math:`\phi`. We construct the centralising element by altering the :math:`\phi`-orbit structure below the orbits of the given *period*.
		
		There are two parameters: a collection of labelled *trees* and a *rearranger* element; in the notation of [BGG11]_ these are elements of :math:`K_{m_i}` and :math:`G_{n, r_i}` respectively.
		.. todo:: doctests
		.. caution:: This is an experimental feature based on [BGG11]_.
		"""
		if period not in self.cycle_type:
			raise ValueError("There are no orbits of period {}.".format(period))
		orbits = self.periodic_orbits[period]
		if not len(orbits) == len(trees):
			raise ValueError("There are {} orbits but {} trees were specified.".format(
				len(orbits), len(trees)
			))
		expected_signature = Signature(self.signature.arity, len(orbits))
		if rearranger.signature != expected_signature:
			raise ValueError("rearranger signature was {}; expected {}.".format(
				rearranger.signature, expected_signature
			))
		
		new_domain = copy(self.domain)
		new_range  = copy(self.range)
		
		for i, orbit in enumerate(orbits):
			tree, shifts = trees[i]
			for j, old_leaf in enumerate(orbit):
				k = new_domain.index(old_leaf)
				
				d_expanded = []
				r_expanded = []
				
				for leaf, shift in zip(tree, shifts):
					tail = leaf[1:]
					d_expanded.append(old_leaf.extend(tail))
					r_expanded.append(self._replacement(
						orbits, i, j, shift, tail, rearranger
					))
				new_domain[k:k+1] = d_expanded
				new_range[k:k+1]  = r_expanded
			
		return Automorphism(new_domain, new_range)
	
	def _replacement(self, orbits, orbit_num, orbit_pos, shift, tail, rearranger):
		period = len(orbits[0])
		new_position = (orbit_pos + shift) % period
		
		relabelling_domain = Generators.standard_basis(rearranger.signature)
		relabelling_range  = Generators(self.signature, (
			orbit[new_position] for orbit in orbits
		))
		relabeller = Homomorphism(relabelling_domain, relabelling_range)
		
		source = Word([orbit_num + 1], rearranger.signature)
		source = source.extend(tail)
		rearrangement = rearranger.image(source)
		target = relabeller.image(rearrangement)
		
		return target


def type_b_triple(power, head, tail):
	return dict(start_tail = tuple(), power=power, end_tail=tail, target=head)

def trapezium_area(x0, y0, x1, y1):
	"""The area of a trapezium whose edges are :math:`x = x_0`, :math:`x = x_1`, :math:`y=x` and the straight line from :math:`(x_0, y_0)` to :math:`(x_1, y_1)`."""
	perp_height = x1 - x0
	avg_width = (abs(y1-x1) + abs(y0-x0)) / 2
	return perp_height * avg_width
