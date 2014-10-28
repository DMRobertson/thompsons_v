r"""As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. 

This module works with automorphisms of :math:`V_{n,r}(\boldsymbol{x})`. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.
 
.. testsetup::
	
	from thompson.word import Word, from_string
	from thompson import word
	from thompson.generators import Generators
	from thompson.automorphism import *
	from thompson.orbits import dump_orbit_types
	from thompson.examples import *
"""

__all__ = ["Automorphism"]

from copy import copy
from io import StringIO
from itertools import product

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
	"""
	def __init__(self, domain, range):
		r"""Creates an automorphism, given the *arity* :math:`n` and *alphabet_size* :math:`r`. Two bases *domain* and *range* are also given. The automorphism maps elements so that the given order of generators is preserved:
		
			.. math:: \text{domain}_{\,i} \mapsto \text{range}_{\,i}
		
		:raises IndexError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis :meth:`isn't actually a basis <thompson.generators.Generators.is_basis>`.
		"""
		#Check to see that the domain and range algebras are the same
		if domain.signature != range.signature:
			raise TypeError('Domain {} and range {} have different signatures.'.format(
			  domain.signature, range.signature))
		
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
		self._inv = {}
		self._qnf_basis = None
		self._qnf_orbit_types = {}
		self.signature = domain.signature
		
		super().__init__(domain, range)
		
		#Setup the inverse map
		for root in Generators.standard_basis(domain.signature):
			self._image_simple_above_domain(root, self.range.signature, self.domain.signature, self._inv)
	
	def _set_image(self, d, r, sig_in, sig_out, cache):
		if cache is not self._map and cache is not self._inv:
			raise ValueError("Incorrect cache provided.")
		super()._set_image(d, r, sig_in, sig_out, cache)
		cache = self._inv if cache is self._map else self._map
		cache[r] = d
	
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
		"""
		#todo examples
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
		
		.. seealso:: Quasi-normal forms are introduced in section 4.2 of the paper. In particular, this method implements Lemma 4.24.1. Higman first described the idea of quasi-normal forms in section 9 of [Hig]_.
		"""
		if self._qnf_basis is not None:
			return self._qnf_basis
		basis = self._seminormal_form_start_point()
		#expand basis until each no element's orbit has finite intersection with X<A>
		i = 0
		while i < len(basis):
			type, image = self._orbit_type(basis[i], basis)
			if type is OrbitType.incomplete:
				basis.expand(i)
			else:
				self._qnf_orbit_types[basis[i]] = type
				i += 1
		self._qnf_basis = basis
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
	
	def _orbit_type(self, y, basis):
		"""Returns the orbit type of *y* with respect to the given *basis*. Also returns a dictionary of computed images, the list (7) from the paper.
		
		>>> #Example 4.5.
		>>> dump_orbit_types(example_4_5, example_4_5.domain)
		x1 a1 a1 a1: Left semi-infinite orbit with characteristic (-1, a1)
		x1 a1 a1 a2: Bi-infinite orbit containing [x1 a1 a2] a1
		x1 a1 a2: Right semi-infinite orbit with characteristic (1, a2)
		x1 a2 a1: Periodic orbit of order 2
		x1 a2 a2: Periodic orbit of order 2
		with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
		>>> dump_orbit_types(example_4_5, basis=example_4_5.domain, words=['x', 'x a1', 'x a2'])
		x: Incomplete orbit
		x a1: Incomplete orbit
		x a2: Incomplete orbit
		with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
		
		>>> #Example 4.11
		>>> dump_orbit_types(example_4_11)
		x1 a1: Left semi-infinite orbit with characteristic (-1, a1)
		x1 a2: Right semi-infinite orbit with characteristic (1, a2)
		with respect to the basis [x1 a1, x1 a2]
		
		>>> #Example 4.12
		>>> basis = example_4_12._seminormal_form_start_point()
		>>> dump_orbit_types(example_4_12, basis)
		x1 a1: Incomplete orbit
		x1 a2: Incomplete orbit
		with respect to the basis [x1 a1, x1 a2]
		>>> basis.expand(0)
		Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2'])
		>>> dump_orbit_types(example_4_12, basis)
		x1 a1 a1: Bi-infinite orbit containing [x1 a2] a2
		x1 a1 a2: Bi-infinite orbit containing [x1 a2] a1
		x1 a2: Incomplete orbit
		with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2]
		>>> basis.expand(2)
		Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
		>>> dump_orbit_types(example_4_12, basis)
		x1 a1 a1: Periodic orbit of order 4
		x1 a1 a2: Periodic orbit of order 4
		x1 a2 a1: Periodic orbit of order 4
		x1 a2 a2: Periodic orbit of order 4
		with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
		
		>>> #Example 4.25
		>>> dump_orbit_types(example_4_25)
		x1 a1: Right semi-infinite orbit with characteristic (1, a1 a1)
		x1 a2 a1: Bi-infinite orbit containing [x1 a1] a1 a2
		x1 a2 a2: Left semi-infinite orbit with characteristic (-1, a1 a1)
		with respect to the basis [x1 a1, x1 a2 a1, x1 a2 a2]
		"""
		
		# print('Computing the orbit type of', y)
		rsi, rpow1, rpow2, rimages = self._test_semi_infinite(y, basis, backward=False)
		# print('RIGHT:', rsi, rpow1, rpow2, [str(x) for x in rimages])
		lsi, lpow1, lpow2, limages  = self._test_semi_infinite(y, basis, backward=True)
		# print('LEFT: ', lsi, lpow1, lpow2, [str(x) for x in limages])
		
		images = {}
		for i, image in enumerate(rimages):
			images[i] = image
		for i, image in enumerate(limages):
			images[-i] = image
		
		if not (lsi or rsi):
			otype = OrbitType.incomplete
			
		elif lsi and not rsi:
			tail = y.test_above(limages[-1])
			assert tail is not None
			assert len(tail) > 0
			assert lpow1 == 0
			characteristic = -lpow2, tail
			otype = OrbitType.semi_infinite(characteristic, backward=True)
		
		elif rsi and not lsi:
			tail = y.test_above(rimages[-1])
			assert tail is not None
			assert len(tail) > 0
			assert rpow1 == 0
			characteristic = (rpow2, tail)
			otype = OrbitType.semi_infinite(characteristic, backward=False)
			
		elif lsi and rsi:
			if limages[-1] == y:
				assert rimages[-1] == y
				assert lpow1 == rpow1 == 0
				assert lpow2 == rpow2
				otype = OrbitType.periodic(lpow2)
			
			else:
				type_b, tail = basis.test_above(rimages[rpow1])
				type_b_data = rpow1, type_b, tail
				otype = OrbitType.complete_infinite(type_b_data)
		
		return otype, images
	
	def _test_semi_infinite(self, y, basis, backward=False):
		r"""Computes the orbit type of *y* under the current automorphism :math:`\psi` with respect to *basis* in the given direction. The process stops when either:
		
			1. an image :math:`y\psi^m` is computed which is not below the basis.
				- infinite: ``False``
				- start: ``0``
				- end: ``m``
				- images: :math:`y, y\psi, \dotsc, y\psi^m`.
			2. two images :math:`y\psi^l, y\psi^m` have been computed, and both start with the same basis element.
				- infinite: ``True``
				- start: ``l``
				- end: ``m``
				- images: :math:`y, y\psi, \dotsc, y\psi^l`.
		
		Returns the tuple *(infinite, start, end, images)*.
		"""
		image = y
		images = [y]
		m = 0
		
		while True:
			m += 1
			image = self.image(image, inverse=backward) #Compute y\phi^{\pm m}
			images.append(image)
			#1. Did we fall out of X<A>?
			if not basis.is_above(image):
				return False, 0, m, images
			
			#2. Is this an infinite orbit?
			#Find the generator *prefix* above the current image.
			prefix, image_tail = basis.test_above(image)
			# print(image, ':', prefix, '|', image_tail)
			for ell, previous in enumerate(images[:-1]):
				if prefix.is_above(previous):
					return True, ell, m, images
	
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
			{..., -1, 0, 1, 2, 3, 4, ...}
		
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
			
		
		:returns: The (possibly empty) :class:`~thompson.orbits.SolutionSet` of all integers :math:`m` for which :math:`u\psi^m = v`. Note that if :math:`u = v` this method returns :math:`m \in \mathbb{Z}`. 
		
		.. seealso:: The implementation is due to lemma 4.24.2 of the paper.
		"""
		#TODO a script which randomly checks examples to verify.
		if u == v:
			return SolutionSet.all_integers
		basis = self.quasinormal_basis()
		# print('QNF basis is', basis)
		
		if not (basis.is_above(u) and basis.is_above(v)):
			depth = max(u.max_depth_to(basis), v.max_depth_to(basis))
			alphas = range(-1, -self.signature.arity - 1, -1)
			solution_set = SolutionSet.all_integers
			
			# For all strings of length *depth* \Gamma made only from alphas:
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
		u_head_type = self._qnf_orbit_types[u_head]
		v_head_type = self._qnf_orbit_types[v_head]
		# print('u = {} | {}\nhead type: {}'.format(
		  # u_head, format(u_tail), u_head_type))
		# print('v = {} | {}\nhead type: {}'.format(
		  # v_head, format(v_tail), v_head_type))
		
		#Is either element periodic?
		if u_head_type.is_type('A') or v_head_type.is_type('A'):
			#Are they both periodic with the same periods?
			if u_head_type != v_head_type:
				# print('Not both periodic with same period.')
				return SolutionSet.empty_set
			
			period = u_head_type.data
			image = u
			# print('Testing the orbit of u.')
			for i in range(1, period):
				image = self.image(image)
				if image == v:
					# print('Found a match.')
					return SolutionSet(i, period)
			# print('No match found.')
			return SolutionSet.empty_set
		
		# print('Neither u nor v is periodic.')
		#Replace them with nicer versions from the same orbits.
		#Track how many positions in the orbit we moved.
		u_shift, u_head, u_tail = self._preprocess(u_head, u_tail)
		v_shift, v_head, v_tail = self._preprocess(v_head, v_tail)
		
		u = u_head.extend(u_tail)
		v = v_head.extend(v_tail)
		# print("Shifted u by {} positions.\nNow u is {} | {}".format(
		  # u_shift, u_head, format(u_tail)))
		# print("Shifted v by {} positions.\nNow v is {} | {}".format(
		  # v_shift, v_head, format(v_tail)))
		
		type, images = self._orbit_type(u, basis)
		for i, image in images.items():
			if image == v:
				return SolutionSet.singleton(u_shift + i - v_shift)
		return SolutionSet.empty_set
	
	#TODO needs a better name
	def _preprocess(self, head, tail):
		r"""Takes a pair :math:`(y, \Gamma)` below the quasi-normal basis :math:`X` and returns a triple :math:`(\widetilde{y}, \widetilde{\Gamma}, k) where
		* The orbit of :math:`\widetilde{y}` is semi-infinite;
		* :math:`\widetilde{y}\widetilde{\Gamma}` is in the same orbit as :math:`y\Gamma`.
		* :math:`\widetilde{\Gamma}` does not start with the characteristic multiplier of :math:`\widetilde{y}`.
		
			>>> example_4_25.quasinormal_basis()
			Generators((2, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
			>>> example_4_25._preprocess(Word('x a2 a2', (2, 1)), from_string('a1 a1 a2'))
			(1, Word('x1 a2 a2', (2, 1)), (-2,))
		"""
		head_type = self._qnf_orbit_types[head]
		if head_type.is_type('C'):
			shift_1, head, prefix = head_type.data
			tail = prefix + tail
		else:
			shift_1 = 0
		
		head_type = self._qnf_orbit_types[head]
		assert head_type.is_type('B')
		characteristic = head_type.data
		shift_2, tail = self._remove_from_start(tail, characteristic)
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
	
	#Testing conjugacy
	def test_conjugate_to(self, other):
		r"""Determines if the current automorphism :math:`\psi` is conjugate to the *other* automorphism :math:`\phi`.
		
		:returns: if it exists, a conjugating element :math:`\rho` such that :math:`\rho^{-1}\psi\rho = \phi`. If no such :math:`\rho` exists, returns ``None``.
		:raises ValueError: if the automorphisms have different arities or alphabet sizes.
		
		.. seealso:: This is an implementation of Algorithm 5.6 in the paper. It depends on Algorithms 5.13 and 5.27. See also the periodic and infinite versions of :meth:`~thompson.factors.PeriodicFactor.test_conjugate_to` and :meth:`~thompson.factors.InfiniteFactor.test_conjugate_to`.
		
			>>> psi = example_5_12_psi; phi = example_5_12_phi
			>>> rho = psi.test_conjugate_to(phi)
			>>> print(rho)
			Automorphism: V(2, 1) -> V(2, 1) specified by 6 generators (after reduction).
			x1 a1 a1 a1 a1 -> x1 a1 a2   
			x1 a1 a1 a1 a2 -> x1 a2 a2   
			x1 a1 a1 a2    -> x1 a1 a1 a1
			x1 a1 a2       -> x1 a2 a1 a1
			x1 a2 a1       -> x1 a1 a1 a2
			x1 a2 a2       -> x1 a2 a1 a2
			>>> rho * phi == psi * rho
			True
		"""
		#todo doctests
		#TODO broken links.
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
			o_i = self.free_factor(o_qnf_i, infinite=True)
			
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
		#todo doctest
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
		periodic = []
		infinite = []
		for gen in basis:
			try:
				type = self._qnf_orbit_types[gen]
			except KeyError:
				type, _ = self._orbit_type(gen, basis)
			if type.is_type("A"):
				periodic.append(gen)
			elif type.is_type("B") or type.is_type("C"):
				infinite.append(gen)
			else:
				raise ValueError('Automorphism was not in semi-normal form with respect to the given basis.')
		periodic = Generators(self.signature, periodic)
		infinite = Generators(self.signature, infinite)
		return periodic, infinite
	
	def free_factor(self, generators, infinite=False):
		r"""This method restricts the current automorphism to the subalgebra generated by the given set :math:`X` of *generators*. This is then transformed into an automorphism of an isomorphic algebra with minimal alphabet size :math:`1 \le s \le n-1`.
		
		.. math::
		
			&G_{n, r}	&\to		&G_{n, |X|}		&\to		&G_{n, s} \\
			&\phi		&\mapsto	&\phi|_X		&\mapsto	&\phi\,'
		
		:raises ValueError: if an empty list of *generators* is provided.
		
		:returns: The transformed automorphism :math:`\phi\, \in G_{n, s}`. Its type is :class:`PeriodicFactor` if *infinite* is False; otherwise its type is :class:`InfiniteFactor`.
		
		.. doctest::
			
			>>> qnb = example_5_3.quasinormal_basis()
			>>> p, i = example_5_3._partition_basis(qnb)
			>>> print(example_5_3.free_factor(p, infinite=False))
			PeriodicFactor: V(2, 1) -> V(2, 1) specified by 2 generators (after reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 a2 a1 ~>    y1 a1 => y1 a2    ~> x1 a1 a2 a2
			x1 a1 a2 a2 ~>    y1 a2 => y1 a1    ~> x1 a1 a2 a1
		
		.. doctest::
			:hide:
			
			>>> # alphabet_size_two example
			>>> qnb = alphabet_size_two.quasinormal_basis()
			>>> p, i = alphabet_size_two._partition_basis(qnb)
			>>> print(alphabet_size_two.free_factor(p, infinite=False))
			PeriodicFactor: V(3, 1) -> V(3, 1) specified by 1 generators (after reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 ~>    y1 => y1    ~> x1 a1
			>>> print(alphabet_size_two.free_factor(i, infinite=True))
			InfiniteFactor: V(3, 1) -> V(3, 1) specified by 5 generators (after reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a2 ~>    y1 a1    => y1 a1 a3    ~> x1 a2 a3
			x1 a3 ~>    y1 a2    => y1 a3       ~> x2      
			x2 a1 ~>    y1 a3 a1 => y1 a2       ~> x1 a3   
			x2 a2 ~>    y1 a3 a2 => y1 a1 a2    ~> x1 a2 a2
			x2 a3 ~>    y1 a3 a3 => y1 a1 a1    ~> x1 a2 a1
		"""
		#TODO more doctests
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
		inverse_relabeller = Homomorphism(copy(generators_relabelled), copy(generators))
		
		#3. Return the factor
		from .factors import get_factor_class
		type = get_factor_class(infinite)
		factor = type(domain, range, inverse_relabeller, inverse_relabeller)
		factor._qnf_basis = copy(generators_relabelled) #See the discussion before ex 5.3
		#TODO pass on the any data about orbit types here
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
		
		domain, range = zip(*sorted(zip(domain, range)))
		domain = Generators(self.signature, domain)
		range = Generators(self.signature, range)
		
		return Automorphism(domain, range)

#TODO? Compose and invert automorphisms. Basically move all the functionality from TreePair over to this class and ditch trree pair.
#TODO method to decide if the automorphism is in (the equivalent of) F, T, or V.
#TODO the named elements A, B, C, X_n of Thompson's V.

