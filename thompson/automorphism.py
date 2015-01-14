r"""As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. 

This module works with automorphisms of :math:`V_{n,r}(\boldsymbol{x})`. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.
 
.. testsetup::
	
	from thompson.word import Word, Signature, from_string, free_monoid_on_alphas
	from thompson import word
	from thompson.generators import Generators
	from thompson.automorphism import *
	from thompson.orbits import print_component_types
	from thompson.examples import *
	from thompson.number_theory import gcd
"""

__all__ = ["Automorphism"]

from builtins import print as builtin_print
from copy import copy
from io import StringIO
from itertools import product, chain

from .word import *
from .generators import *
from .homomorphism import Homomorphism
from .orbits import *

def modulo_non_zero(x, n):
	r"""Returns the unique integer :math:`s` such that :math:`1 \le s \le n` and :math:`s \equiv x \pmod{n}`.
	
		>>> [modulo_non_zero(x, 10) for x in range(20)]
		[10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
	"""
	x %= n
	if x == 0:
		return n
	return x

class Automorphism(Homomorphism):
	"""
	:ivar signature: The :class:`~thompson.word.Signature` shared by the generating sets domain and range.
	:ivar pond_banks: A list of tuples :math:`(\ell, k, r)` such that :math:`(\ell, r)` are banks of a pond with :math:`\ell\phi^k = r`.
	"""
	#todo docstring for reduce
	def __init__(self, domain, range, reduce=True):
		r"""Creates an automorphism, given the *arity* :math:`n` and *alphabet_size* :math:`r`. Two bases *domain* and *range* are also given. The automorphism maps elements so that the given order of generators is preserved:
		
			.. math:: \text{domain}_{\,i} \mapsto \text{range}_{\,i}
		
		:raises TypeError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis :meth:`isn't actually a basis <thompson.generators.Generators.is_basis>`.
		"""
		#Check to see that the domain and range algebras are the same
		if domain.signature != range.signature:
			raise TypeError('Domain {} and range {} have different signatures.'.format(
			  domain.signature, range.signature))
		
		super()._expand(domain, range)
		
		#Check to see that range is a basis for the range algebra
		i, j = range.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{}. Missing elements are {}.".format(
			  range.signature, [format(x) for x in missing]))
		
		#Cache attributes
		self.pond_banks = None
		self._inv = {}
		self._qnf_basis = None
		self.signature = domain.signature
		self._characteristics = None
		
		super().__init__(domain, range, reduce)
		
		#Setup the inverse map
		for root in Generators.standard_basis(domain.signature):
			self._image_simple_above_domain(root, self.range.signature, self.domain.signature, self._inv)
	
	def _set_image(self, d, r, sig_in, sig_out, cache):
		if cache is not self._map and cache is not self._inv:
			raise ValueError("Incorrect cache provided.")
		super()._set_image(d, r, sig_in, sig_out, cache)
		cache = self._inv if cache is self._map else self._map
		cache[r] = d
	
	def __invert__(self):
		"""We overload python's unary negation operator ``~`` as shorthand for inversion. (In Python, ``~`` is normally used for bitwise negation.) We can also call a method explicitily: ``phi.inverse()`` is exactly the same as ``~phi``.
		
			>>> phi = random_automorphism(20, signature = Signature(2, 1))
			>>> phi * ~phi == ~phi * phi
			True
			>>> print(phi * ~phi)
			Automorphism: V(2, 1) -> V(2, 1) specified by 1 generators (after expansion and reduction).
			x1 -> x1
			
		"""
		inv = copy(self)
		inv.domain, inv.range = Generators.sort_mapping_pair(self.range, self.domain)
		inv._map, inv._inv = self._inv, self._map
		inv._qnf_basis = None
		return inv
		
	inverse = __invert__
	
	@classmethod
	def identity(cls, signature):
		"""Creates a new automorphism which is the identity map on the algebra with the given *signature*.
		
			>>> print(Automorphism.identity((3, 2)))
			Automorphism: V(3, 2) -> V(3, 2) specified by 2 generators (after expansion and reduction).
			x1 -> x1
			x2 -> x2
		"""
		d = Generators.standard_basis(signature)
		r = Generators.standard_basis(signature)
		return cls(d, r)
	
	def is_identity(self):
		"""Returns True if this automorphism is the identity map on the algebra with the given *signature*. Otherwise returns False.
		
			>>> aut = Automorphism.identity(random_signature())
			>>> aut.is_identity()
			True
			>>> example_4_25.is_identity()
			False
		"""
		basis = Generators.standard_basis(self.signature)
		return self.image_of_set(basis) == basis
	
	#Computing images
	def image(self, key, inverse=False):
		"""If *inverse* is True, the inverse of the current automorphism is used to map *key* instead. Otherwise this method delegates to :meth:`Homomorphism.image <thompson.homomorphism.Homomorphism.image>`.
		
		Examples of finding inverse images:
		
			>>> print(example_4_25.image('x1 a2 a2', inverse=True))
			x1 a2 a2 a1 a1
			>>> print(example_4_25.image('x1 a1 a1 a2 a2 a1', inverse=True))
			x1 a2 a1 a2 a1
			>>> print(example_4_25.image('x a2', inverse=True))
			x1 a2 a2 a2 x1 a2 a2 a1 a1 L
			>>> print(example_4_25.image('x a2 a2 x a1 a2 L', inverse=True))
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
			>>> print(example_5_3.image_of_set(basis, inverse=True))
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
			
			>>> print(example_4_25.repeated_image('x1 a1', 10))
			x1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1
			>>> print(example_4_25.repeated_image('x1 a1 a1 a1 a1 a1 a1 a1', -3))
			x1 a1
			>>> print(arity_four.repeated_image('x1 a4 a4 a2', 4))
			x1 a3 a3 a2
			>>> print(arity_four.repeated_image('x1 a3 a3 a2', -4))
			x1 a4 a4 a2
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
	
	#Generating the quasinormal basis
	#Some ideas for making this process faster:
	#	Cache the component types (wrt QNB) for the QNB elements
	#	Put some part of the core of type B orbits into confirmed? See /theory/Characteristic components and the QNB.tex
	#	Cache the powers (-ell, m) which determine the core part of a component
	def quasinormal_basis(self):
		r"""We say that :math:`\phi` is *in semi-normal form* with respect to the basis :math:`X` if no element of :math:`X` lies in an incomplete :math:`X`-component of a :math:`\phi` orbit. See the :mod:`~thompson.orbits` module for more details.
		
		There is a minimal such basis, :math:`X_\phi` say, and we say that :math:`\phi` is *in quasi-normal form* with respect to :math:`X_\phi`. This method determines and returns the basis :math:`X_\phi` where :math:`\phi` denotes the current automorphism. The result is cached so that further calls to this method perform no additional computation. 
		
			>>> example_4_5.quasinormal_basis()
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> alphabet_size_two.quasinormal_basis()
			Generators((3, 2), ['x1 a1', 'x1 a2', 'x1 a3', 'x2'])
			>>> example_5_12_phi.quasinormal_basis()
			Generators((2, 1), ['x1 a1', 'x1 a2'])
		
		:rtype: a :class:`~thompson.generators.Generators` instance.
		
		.. note:: This method additionally looks for ponds and caches the data describing any ponds it finds.
		
		.. seealso:: Quasi-normal forms are introduced in section 4.2 of the paper. In particular, this method implements Lemma 4.24.1. Higman first described the idea of quasi-normal forms in section 9 of [Hig]_.
		"""
		#0. Have we already computed the basis?
		if self._qnf_basis is not None:
			return self._qnf_basis
		
		#1. Expand the starting basis until each no element's belongs to a finite X-component.
		confirmed = set()
		basis = self._seminormal_form_start_point()
		basis.cache = set(basis)
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
			else:
				if ctype.is_type_A():
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
				
				elif ctype.is_type_B():
					confirmed.add(basis[i])
				i = (i + 1) % len(basis)
		
		self._qnf_basis = basis
		
		#2. Look for and remember the details of any pond orbits.
		terminal, initial = self.semi_infinite_end_points()
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
		
		self.pond_banks = ponds
		return basis
	
	def _seminormal_form_start_point(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \cup W`. See remark 4.10 and example 4.25. This basis that this method produces is the smallest possible which *might* be semi-normal.
		
			>>> example_4_5._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> example_4_11._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1', 'x1 a2'])
			>>> example_4_12._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1', 'x1 a2'])
			>>> example_4_25._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
			>>> cyclic_order_six._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2 a1', 'x1 a1 a2 a2', 'x1 a2'])
		"""
		basis = Generators.standard_basis(self.signature)
		i = 0
		while i < len(basis):
			b = basis[i]
			if b in self.domain or b in self.range:
				i += 1
			else:
				basis.expand(i)
		return basis
	
	def orbit_type(self, y, basis):
		"""Returns the orbit type of *y* with respect to the given *basis*. Also returns a dictionary of computed images, the list (7) from the paper.
		
			>>> #Example 4.5.
			>>> print_component_types(example_4_5, example_4_5.domain)
			x1 a1 a1 a1: Left semi-infinite component with characteristic (-1, a1)
			x1 a1 a1 a2: Bi-infinite component
			x1 a1 a2: Right semi-infinite component with characteristic (1, a2)
			x1 a2 a1: Periodic component of order 2
			x1 a2 a2: Periodic component of order 2
			with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> print_component_types(example_4_5, basis=example_4_5.domain, words=['x', 'x a1', 'x a2'])
			x1: Incomplete component
			x1 a1: Incomplete component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			
			>>> #Example 4.11
			>>> print_component_types(example_4_11)
			x1 a1: Left semi-infinite component with characteristic (-1, a1)
			x1 a2: Right semi-infinite component with characteristic (1, a2)
			with respect to the basis [x1 a1, x1 a2]
			
			>>> #Example 4.12
			>>> basis = example_4_12._seminormal_form_start_point()
			>>> print_component_types(example_4_12, basis)
			x1 a1: Incomplete component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1, x1 a2]
			>>> basis.expand(0)
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2'])
			>>> print_component_types(example_4_12, basis)
			x1 a1 a1: Bi-infinite component
			x1 a1 a2: Bi-infinite component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2]
			>>> basis.expand(2)
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> print_component_types(example_4_12, basis)
			x1 a1 a1: Periodic component of order 4
			x1 a1 a2: Periodic component of order 4
			x1 a2 a1: Periodic component of order 4
			x1 a2 a2: Periodic component of order 4
			with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			
			>>> #Example 4.25
			>>> print_component_types(example_4_25)
			x1 a1: Right semi-infinite component with characteristic (1, a1 a1)
			x1 a2 a1: Bi-infinite component
			x1 a2 a2: Left semi-infinite component with characteristic (-1, a1 a1)
			with respect to the basis [x1 a1, x1 a2 a1, x1 a2 a2]
		
		.. seealso:: Lemmas 4.14, 4.24
		"""
		images = {}
		type_b_data = None
		
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
				characteristic = -lpow2, tail
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
				characteristic = rpow2, tail
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
			
			2. For some :math:`0 \le l < m`, :math:`y\psi^l` and :math:`y\psi^m` start with the same basis element.
			
				- infinite: ``True``
				- start: ``l``
				- end: ``m``
				- images: :math:`y, y\psi, \dotsc, y\psi^{m}`.
				
			.. note::
				
				The word :math:`y\psi^m` is not strictly in the core part of the orbit of Lemma 4.24. We return this as part of *images* so that we can compute the characteristic multiplier in :meth:`orbit_type`.
		
		:returns: the tuple *(infinite, start, end, images)*.
		
		.. seealso:: Lemma 4.24
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
					if prefix == previous[:len(prefix)]:
						images.append(image)
						return True, ell, m, images
			images.append(image)
	
	def semi_infinite_end_points(self):
		r"""Returns the list of terminal :class:`Words <thompson.word.Word>` in left semi-infinite components and the list of initial words in right semi-infinite components. This is all computed with respect the current automorphism's :meth:`quasinormal_basis`. These are the sets :math:`X\langle A\rangle \setminus Y\langle A\rangle` and :math:`X\langle A\rangle \setminus Z\langle A\rangle`.
		
			>>> print(*example_4_5.semi_infinite_end_points())
			[x1 a1 a1] [x1 a1 a2]
			>>> print(*example_4_11.semi_infinite_end_points())
			[x1 a1] [x1 a2]
			>>> print(*example_4_12.semi_infinite_end_points())
			[] []
			>>> print(*example_4_25.semi_infinite_end_points())
			[x1 a2 a2, x1 a2 a2 a1] [x1 a1, x1 a1 a1]
		
		:rtype: A pair of :class:`Generators <thompson.generators.Generators>`.
		
		.. seealso:: The discussion before lemma 4.6.
		"""
		basis = self.quasinormal_basis()                  #X
		min_expansion = basis.minimal_expansion_for(self) #Y
		img_expansion = self.image_of_set(min_expansion)  #Z
		terminal = basis.descendants_above(min_expansion) #X<A> \ Y<A>
		initial  = basis.descendants_above(img_expansion) #X<A> \ Z<A>
		return terminal, initial
	
	def _descend_to_complete_infinite(self, endpt):
		"""A really crude use of the lemma described in AJD's email.
		
		**Lemma.** Let :math:`b` be the bank of a pond. There is a :math:`\Gamma \in A^*` such that :math:`b\Gamma` belongs to a doubly-infinite :math:`\phi`-orbit that meets :math:`X`.
		"""
		basis = self.quasinormal_basis()
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
	
	def characteristics(self, print=False):
		"""This method computes the set of characteristics :math:`(m, \Gamma)` of all characteristic words with respect to the :meth:`quasinormal_basis`.
		
		:param bool print: if True, prints the characteristics in a nicely formatted way and does **not** return anything.
		
		:rtype: a :class:`set`.
		
		.. doctest::
			
			>>> from pprint import pprint
			>>> example_4_25.characteristics(print=True)
			(-1, a1 a1)
			(1, a1 a1)
			>>> psi, phi = random_conjugate_pair()
			>>> #Lemma 5.16
			>>> psi.characteristics() == phi.characteristics()
			True
			>>> if not psi.characteristics() == phi.characteristics():
			... 	print(psi.characteristics())
			... 	print(phi.characteristics())
			... 	psi.save_to_file('psi.aut')
			... 	phi.save_to_file('phi.aut')
			>>> (example_6_2 * example_6_2).characteristics(print=True)
			(-1, a1)
			(1, a2 a2)
		
		.. doctest::
			:hide:
			
			>>> #Lemma 6.1
			>>> from random import randint
			>>> psi = random_conjugate_infinite_factors()[0]
			>>> original_chars = psi.characteristics()
			>>> power = psi
			>>> a = randint(2, 6)
			>>> for _ in range(a - 1): 
			... 	power *= psi
			>>> chars = set()
			>>> for mult, char in original_chars:
			... 	d = gcd(mult, a)
			... 	q = abs(a // d)
			... 	chars.add((mult//d, char*q))
			>>> chars == power.characteristics()
			True
		
		.. seealso:: Defintion 5.14.
		"""
		if self._characteristics is None:
			basis = self.quasinormal_basis()
			self._characteristics = set()
			initial, terminal = self.semi_infinite_end_points()
			for endpt in chain(initial, terminal):
				ctype, _, _ = self.orbit_type(endpt, basis)
				if ctype.is_type_B():
					self._characteristics.add(ctype.characteristic)
		if print:
			for power, mult in sorted(self._characteristics):
				builtin_print('({}, {})'.format(power, format(mult)))
		else:
			return self._characteristics
	
	#Orbit sharing test
	def share_orbit(self, u, v):
		r"""Determines if :math:`u` and :math:`v` are in the same orbit of the current automorphism :math:`\psi`. Specifically, does there exist an integer :math:`m` such that :math:`u\psi^m = v`?
		
		.. doctest::
			:options: -ELLIPSIS
			
			>>> u  = Word('x1 a2 a3 a1 a2', (3, 2))
			>>> v1 = Word('x1 a1 a2 a2 a3 a1', (3, 2))
			>>> v2 = Word('x2 a3 a2', (3, 2))
			>>> print(alphabet_size_two.share_orbit(u, v1))
			{}
			>>> print(alphabet_size_two.share_orbit(u, v2))
			{-2}
			>>> print(alphabet_size_two.share_orbit(u, u))
			{0}
		
		.. doctest::
			:hide:
			:options: -ELLIPSIS
			
			>>> u  = Word('x a2 a2 a1 a1 a2', (2, 1))
			>>> v1 = Word('x a1 a2', (2, 1))
			>>> v2 = Word('x a1 a1 a2', (2, 1))
			>>> print(example_4_25.share_orbit(u, v1))
			{}
			>>> print(example_4_25.share_orbit(u, v2))
			{3}
			>>> u  = Word('x a2 a2 x a1 a2 L x a2 a1 L x a1 a1 a1 a2 L', (2, 1))
			>>> vs = [
			... 	Word('x a2 a2 a1 a1 a1 a1 a1', (2, 1)),
			... 	Word('x a2 a2 a2', (2, 1)),
			... 	Word('x a2 a2 x a1 a2 L', (2, 1)),
			... 	Word('x a1 a1 x a1 a2 x a2 a2 L L', (2, 1))]
			... 
			>>> for v in vs: print(example_4_25.share_orbit(u, v))
			{-4}
			{}
			{-1}
			{}
			>>> u1 = Word('x a2 a3 a1', (4, 1))
			>>> v1 = Word('x a3 a3 a3 a3 a3 a3 a3 a3 a3 a2 a3 a1', (4, 1))
			>>> v2 = Word('x a1 a2 a3 a4', (4, 1))
			>>> u2 = Word('x a3 a4 a1', (4, 1))
			>>> v3 = Word('x a4 a1 a1', (4, 1))
			>>> v4 = Word('x a4 a3 a2 a1', (4, 1))
			>>> print(arity_four.share_orbit(u1, v1))
			{9}
			>>> print(arity_four.share_orbit(u1, v2))
			{}
			>>> print(arity_four.share_orbit(u2, v3))
			{-1}
			>>> print(arity_four.share_orbit(u2, v4))
			{}
			>>> u = Word("x a1 a2 a1 a2 a1", (2, 1))
			>>> v = Word("x a2 a2 a2 a1", (2, 1))
			>>> print(cyclic_order_six.share_orbit(u, v))
			{..., -1, 2, 5, 8, 11, 14, ...}
			>>> u = Word("x a1 a1 x a1 a2 a1 x a1 a2 a2 a1 L L", (2, 1))
			>>> v1 = Word("x a1 a1 x a2 a2 x a2 a1 L L", (2, 1))
			>>> v2 = Word("x a1 a1 x a1 a2 a1 x a1 a1 L L", (2, 1))
			>>> print(cyclic_order_six.share_orbit(u, v1))
			{..., -1, 5, 11, 17, 23, 29, ...}
			>>> print(cyclic_order_six.share_orbit(u, v2))
			{}
			>>> #Two sides of a pond
			>>> u = Word('x a1 a1 a1 a1 a1 a1 a1 a2', (2, 1))
			>>> v = Word('x a2 a2 a1 a1 a2', (2, 1))
			>>> print(first_pond_example_phi.share_orbit(u, v))
			{4}
		
		:returns: The (possibly empty) :class:`~thompson.orbits.SolutionSet` of all integers :math:`m` for which :math:`u\psi^m = v`. Note that if :math:`u = v` this method returns :math:`\mathbb{Z}`. 
		
		.. seealso:: The implementation is due to lemma 4.24.2 of the paper.
		"""
		#TODO a script which randomly checks examples to verify.
		#TODO pond examples
		orig_u, orig_v = u, v
		basis = self.quasinormal_basis()
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
		
			>>> basis = example_4_25.quasinormal_basis()
			>>> basis
			Generators((2, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
			>>> head = Word('x a2 a2', (2, 1))
			>>> _, _, type_b_data = example_4_25.orbit_type(head, basis)
			>>> example_4_25._type_b_descendant(head, from_string('a1 a1 a2'), type_b_data)
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
		
		head_type, _, _ = self.orbit_type(head, self.quasinormal_basis())
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
		basis = self.quasinormal_basis()
		
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
	
	#Testing conjugacy
	def test_conjugate_to(self, other):
		r"""Determines if the current automorphism :math:`\psi` is conjugate to the *other* automorphism :math:`\phi`.
		
		:returns: if it exists, a conjugating element :math:`\rho` such that :math:`\rho^{-1}\psi\rho = \phi`. If no such :math:`\rho` exists, returns ``None``.
		:raises ValueError: if the automorphisms have different arities or alphabet sizes.
		
		.. doctest::
			
			>>> psi, phi = random_conjugate_pair()
			>>> rho = psi.test_conjugate_to(phi)
			>>> rho is not None
			True
			>>> psi * rho == rho * phi
			True
			>>> psi, phi = first_pond_example_psi, first_pond_example_phi
			>>> rho = psi.test_conjugate_to(phi)
			>>> rho is not None
			True
			>>> psi * rho == rho * phi
			True
		
		..doctest::
			:hide:
			
			>>> rho = nathan_pond_example.test_conjugate_to(nathan1_example)
			>>> rho is None
			True
			>>> rho = nathan2_example.test_conjugate_to(nathan3_example)
			>>> rho is not None
			True
			>>> nathan2_example * rho == rho * nathan3_example
			True
			>>> rho = nathan4_example
			>>> nathan2_example * rho == rho * nathan3_example
			True
			
		
		.. seealso:: This is an implementation of Algorithm 5.6 in the paper. It depends on Algorithms 5.13 and 5.27; these are the :meth:`periodic <thompson.factors.PeriodicFactor.test_conjugate_to>` and :meth:`infinite <thompson.factors.InfiniteFactor.test_conjugate_to>` tests for conjugacy.
		"""
		#TODO Doctest: try assembling a conjugator from factors
		#0. Check that both automorphisms belong to the same G_{n,r}.
		if self.signature != other.signature:
			raise ValueError('Automorphism\'s signatures {} and {} do not match.'.format(
			  self.signature, other.signature))
		
		#1. Before going any further, check that the number of periodic and infinite elements are compatible.
		result = self._check_parition_sizes(other)
		if result is None:
			return None
		pure_periodic, pure_infinite, s_qnf_p, s_qnf_i, o_qnf_p, o_qnf_i = result
		
		#4. If necessary, test the periodic factors.
		if pure_infinite:
			rho_p = None
		else:
			#2, 3. Extract the periodic factor.
			s_p = self.free_factor(s_qnf_p, infinite=False)
			o_p = other.free_factor(o_qnf_p, infinite=False)
			rho_p = s_p.test_conjugate_to(o_p)
			if rho_p is None:
				return None
		
		#Step 5. If necessary, test the infinite factors.
		if pure_periodic:
			rho_i = None
		else:
			#2, 3. Extract the infinite factor.
			s_i = self.free_factor(s_qnf_i, infinite=True)
			o_i = other.free_factor(o_qnf_i, infinite=True)
			rho_i = s_i.test_conjugate_to(o_i)
			if rho_i is None:
				return None
		
		#Step 6. If we've got this far, we have a periodic/infinite mixture.
		#Combine the rewritten conjugators into one conjugator in G_n,r
		return self._combine_factors(rho_p, rho_i)
	
	def _check_parition_sizes(self, other):
		r"""Checks the sizes of a partition of two quasinormal bases to see if they might possibly describe conjugate automorphisms.
		
		:returns: None if we discover conjugacy is impossible. Otherwise, returns a pair of booleans *(pure_periodic, pure_infinite)*. Note that even if such a pair is returned, we do **not** know for certain that conjugacy is possible.
		"""
		s_qnf = self.quasinormal_basis()
		o_qnf = other.quasinormal_basis()
		#Check that we have at at least one element in each basis.
		#This eliminates the possibility that either s or o is both pure periodic and pure infinite.
		assert (len(s_qnf) != 0 and len(o_qnf) != 0), "One of the QNF bases consists of 0 elements."
		
		s_qnf_p, s_qnf_i = self._partition_basis(s_qnf)
		o_qnf_p, o_qnf_i = other._partition_basis(o_qnf)
		
		s_pure_periodic = len(s_qnf_i) == 0
		s_pure_infinite = len(s_qnf_p) == 0
		o_pure_periodic = len(o_qnf_i) == 0
		o_pure_infinite = len(o_qnf_p) == 0
		
		#Check that s is periodic iff o is periodic and the same for infinite.
		if not (s_pure_periodic == o_pure_periodic
		  and   s_pure_infinite == o_pure_infinite):
			return None
		
		#Check that the lengths match up modulo n-1.
		modulus = self.signature.arity - 1
		if not (len(s_qnf_p) % modulus == len(o_qnf_p) % modulus
		  and   len(s_qnf_i) % modulus == len(o_qnf_i) % modulus):
			return None
		return s_pure_periodic, s_pure_infinite, s_qnf_p, s_qnf_i, o_qnf_p, o_qnf_i
	
	def _partition_basis(self, basis):
		r"""Let the current automorphism be in semi-normal form with respect to the given *basis*. This method places the elements of *basis* into two lists *(periodic, infinite)* depending on their orbit type. 
		
		:returns: the pair *(periodic, infinite)*. Both entries are sets of :class:`~thompson.generators.Generators`.
		
		.. doctest::
			
			>>> basis = example_4_25.quasinormal_basis()
			>>> print(*example_4_25._partition_basis(basis), sep='\n')
			[]
			[x1 a1, x1 a2 a1, x1 a2 a2]
			>>> basis = example_4_5.quasinormal_basis()
			>>> print(*example_4_5._partition_basis(basis), sep='\n')
			[x1 a2 a1, x1 a2 a2]
			[x1 a1 a1, x1 a1 a2]
			>>> basis = example_5_3.quasinormal_basis()
			>>> print(*example_5_3._partition_basis(basis), sep='\n')
			[x1 a1 a2 a1, x1 a1 a2 a2]
			[x1 a1 a1 a1, x1 a1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> basis = example_5_3.quasinormal_basis()
			>>> basis = basis.minimal_expansion_for(example_5_3)
			>>> print(*example_5_3._partition_basis(basis), sep='\n')
			[x1 a1 a2 a1, x1 a1 a2 a2]
			[x1 a1 a1 a1 a1, x1 a1 a1 a1 a2, x1 a1 a1 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2]
		"""
		periodic = Generators(self.signature)
		infinite = Generators(self.signature)
		for gen in basis:
			type, _, _ = self.orbit_type(gen, basis)
			if type.is_incomplete():
				raise ValueError('Automorphism is not in semi-normal form with respect to the given basis.')
			elif type.is_type_A():
				periodic.append(gen)
			else:
				infinite.append(gen)
		return periodic, infinite
	
	def free_factor(self, generators, infinite=False):
		r"""This method restricts the current automorphism to the subalgebra generated by the given set :math:`X` of *generators*. This is then transformed into an automorphism of an isomorphic algebra with minimal alphabet size :math:`1 \le s \le n-1`.
		
		.. math::
		
			&G_{n, r}	&\to		&G_{n, |X|}		&\to		&G_{n, s} \\
			&\phi		&\mapsto	&\phi|_{\langle X\rangle}		&\mapsto	&\phi\,'
		
		:returns: The transformed automorphism :math:`\phi\, \in G_{n, s}`. Its type is :class:`PeriodicFactor` if *infinite* is False; otherwise its type is :class:`InfiniteFactor`.
		:raises ValueError: if an empty list of *generators* is provided.
		
		.. doctest::
			:hide:
			
			>>> qnb = example_5_3.quasinormal_basis()
			>>> p, i = example_5_3._partition_basis(qnb)
			>>> print(example_5_3.free_factor(p, infinite=False))
			PeriodicFactor: V(2, 1) -> V(2, 1) specified by 2 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 a2 a1 ~>    y1 a1 => y1 a2    ~> x1 a1 a2 a2
			x1 a1 a2 a2 ~>    y1 a2 => y1 a1    ~> x1 a1 a2 a1
		
		.. doctest::
			
			>>> # alphabet_size_two example
			>>> qnb = alphabet_size_two.quasinormal_basis()
			>>> p, i = alphabet_size_two._partition_basis(qnb)
			>>> print(alphabet_size_two.free_factor(p, infinite=False))
			PeriodicFactor: V(3, 1) -> V(3, 1) specified by 1 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 ~>    y1 => y1    ~> x1 a1
			>>> print(alphabet_size_two.free_factor(i, infinite=True))
			InfiniteFactor: V(3, 1) -> V(3, 1) specified by 5 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a2 ~>    y1 a1    => y1 a1 a3    ~> x1 a2 a3
			x1 a3 ~>    y1 a2    => y1 a3       ~> x2      
			x2 a1 ~>    y1 a3 a1 => y1 a2       ~> x1 a3   
			x2 a2 ~>    y1 a3 a2 => y1 a1 a2    ~> x1 a2 a2
			x2 a3 ~>    y1 a3 a3 => y1 a1 a1    ~> x1 a2 a1
		"""
		if len(generators) == 0:
			raise ValueError('Must provide at least one generator.')
		#1. Decide how to relabel *generators* as elements of V_n,s
		modulus = self.signature.arity - 1
		alphabet_size = modulo_non_zero(len(generators), modulus)
		generators_relabelled = Generators.standard_basis((self.signature.arity, alphabet_size))
		generators_relabelled.expand_to_size(len(generators))
		
		#2. Relabel the domain and range
		domain = generators.minimal_expansion_for(self)
		range = self.image_of_set(domain)
		
		rewrite_set(domain, generators, generators_relabelled)
		rewrite_set(range,  generators, generators_relabelled)
		#Make sure we can undo the relabelling
		inverse_relabeller = Homomorphism(copy(generators_relabelled), copy(generators), reduce=False)
		
		#3. Return the factor
		from .factors import get_factor_class
		type = get_factor_class(infinite)
		factor = type(domain, range, inverse_relabeller, inverse_relabeller, reduce=False)
		#TODO We know a priori what the quasi-normal basis is. We could pass this info on to the factor?
		#TODO Would also be able to relabel type b data for bi-infinite orbit types
		return factor
	
	def _combine_factors(self, periodic, infinite):
		#doctest and docstring
		if periodic is None:
			p_domain = p_range = Generators(self.signature)
		else:
			p_domain, p_range = periodic.relabel()
		
		if infinite is None:
			i_domain = i_range = Generators(self.signature)
		else:
			i_domain, i_range = infinite.relabel()
		
		assert p_domain.signature == p_range.signature == i_domain.signature == i_range.signature == self.signature
		domain = p_domain + i_domain
		range  = p_range + i_range
		
		domain = Generators(self.signature, domain)
		range = Generators(self.signature, range)
		
		return Automorphism(domain, range)
	
	#For convenience
	def free_factors(self):
		"""A convenience method that produces the periodic and infinite components of an automorphism."""
		basis = self.quasinormal_basis()
		p, i = self._partition_basis(basis)
		p_factor = handle_trivial_factors(self, p, False)
		i_factor = handle_trivial_factors(self, i, True)
		return p_factor, i_factor
	
	def dump_QNB(self):
		basis = self.quasinormal_basis()
		for gen in basis:
			ctype, _, _ = self.orbit_type(gen, basis)
			print(gen, ctype)
	
	#Power conjugacy
	def test_power_conjugate(self, other):
		r"""Determines if some power of the current automorphism :math:`\psi` is conjugate to some power of the *other* automorphism :math:`\phi`.
		
		:returns: if it exists, a triple :math:`(a, b, \rho)` such that :math:`\rho^{-1}\psi^a\rho = \phi^b`. If no such triple exists, returns ``None``.
		:raises ValueError: if the automorphisms have different arities or alphabet sizes.
		
		.. seealso:: This is an implementation of Algorithm 6.12 in the paper. It depends on Algorithms 5.6 (the :meth:`conjugacy test <test_conjugate_to>`) and 6.10 (the :meth:`infinite power conjugate test <thompson.factors.InfiniteFactor.find_power_conjugators>.`)
		"""
		#0. Check that both automorphisms belong to the same G_{n,r}.
		if self.signature != other.signature:
			raise ValueError('Automorphism\'s signatures {} and {} do not match.'.format(
			  self.signature, other.signature))
		
		#1. Before going any further, check that the number of periodic and infinite elements are compatible.
		result = self._check_parition_sizes(other)
		if result is None:
			return None
		pure_periodic, pure_infinite, s_qnf_p, s_qnf_i, o_qnf_p, o_qnf_i = result
		
		#2. Construct PC_p, the set of conjugators that act on the periodic factor
		if not pure_infinite:
			#Extract the periodic factors
			s_p = self.free_factor(s_qnf_p, infinite=False)
			o_p = other.free_factor(o_qnf_p, infinite=False)
			
			#Preprare to brute force search
			bounds = s_p.power_conjugacy_bounds(o_p)
			periodic_conjugators = []
			for data in s_p._test_power_conjugate_upto(o_p, *bounds, inverses=False):
				if pure_periodic:
					a, b, rho_p = data
					return a, b, self._combine_factors(rho_p, None)
				periodic_conjugators.append(data)
			if len(periodic_conjugators) == 0:
				return None
		
		#Step 4. Construct PC_I, the set of conjugators that act on the infinite factor
		if not pure_periodic:
			#Extract the periodic factors
			s_i = self.free_factor(s_qnf_i, infinite=True)
			o_i = other.free_factor(o_qnf_i, infinite=True)
			
			#Preprare to brute force search
			bounds = s_i.power_conjugacy_bounds(o_i)
			infinite_conjugators = []
			for data in s_i._test_power_conjugate_upto(o_i, *bounds, inverses=True):
				if pure_infinite:
					a, b, rho_i = data
					return a, b, self._combine_factors(None, rho_i)
				infinite_conjugators.append(data)
			if len(infinite_conjugators) == 0:
				return None
		
		#5. Try to recombine.
		for alpha, beta, rho_i in infinite_conjugators:
			for c, d, rho_p in periodic_conjugators:
				solns = solve_linear_congruence(a, c, s_p.order) & solve_linear_congruence(b, d, o_p.order)
				if not solns.is_empty():
					soln = solns.base % lcm(s_p.order, o_p.order)
					rho = self._combine_factors(rho_p, rho_i)
					return alpha*soln, beta*soln, rho
		
		#6. If we've got this far, we're out of luck.
		return None
		
	def _test_power_conjugate_upto(self, other, sbound, obound, inverses=False):
		r"""In both the periodic and infinite cases, we establish bounds on the powers :math:`a, b` for conjugacy; the rest is brute force.  This method tests to see if :math:`\psi^a` is conjugate to :math:`\phi^b` within the supplied bounds. Should it find a conjugator :math:`\rho`, this method yields a triple :math:`(a, b, \rho)`.
		
		Let :math:`\hat a, \hat b` denote *sbound* and *obound* respectively. If *inverses* is False, then we search over the ranges :math:`1 \le a \le \hat a` and :math:`1 \le b \le \hat b`. If *inverses* is True, we search over the (four times larger) range `1 \le |a| \le \hat a` and :math:`1 \leq |b| \le \hat b`
		"""
		if sbound < 0:
			raise ValueError('sbound parameter should be at least 0 (received {}).'.format(sbound))
		if obound < 0:
			raise ValueError('obound parameter should be at least 0 (received {}).'.format(obound))
		
		if sbound == 0 or obound == 0:
			return
		
		#Assuming that computation is more expensive than storage space.
		s_powers = dict()
		s_powers[1] = self
		if inverses:
			s_powers[-1] = ~self
		
		o_powers = dict()
		o_powers[1] = other
		if inverses:
			o_powers[-1] = ~other
		
		for i in range(2, sbound + 1):
			s_powers[i] = s_powers[i-1] * self
			if inverse:
				s_powers[-i] = s_powers[1-i] * s_powers[-1]
		
		for i in range(2, obound + 1):
			o_powers[i] = o_powers[i-1] * other
			if inverse:
				o_powers[-i] = o_powers[1-i] * o_powers[-1]
		
		for a, spow in s_powers.items():
			for b, opow in o_powers.items():
				rho = spow.test_conjugate_to(opow)
				if rho is not None:
					yield a, b, rho

def handle_trivial_factors(aut, gens, infinite=False):
	try:
		factor = aut.free_factor(gens, infinite)
	except ValueError as e:
		if e.args[0] == 'Must provide at least one generator.':
			return None
		else:
			raise e
	return factor

def type_b_triple(power, head, tail):
	return dict(start_tail = tuple(), power=power, end_tail=tail, target=head)