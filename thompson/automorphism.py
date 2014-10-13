r"""The algebra of :mod:`words <thompson.word>` has its own structure, and just like groups, rings etc. there exist maps (homomorphisms) which preserve this structure. In our specific case, a homomorphism :math:`\phi: V_{n,r} \to V_{n,s}` is a function satisfying

.. math:: w \alpha_i \phi = w \phi \alpha_i \qquad \text{and} \qquad
          w_1 \dots w_n \lambda \phi = w_1 \phi \dots w_n \phi \lambda.

As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. 

This module works with automorphisms of :math:`V_{n,r}(\boldsymbol{x})`. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.
 
.. testsetup::
	
	from thompson.automorphism import *
	from thompson.word import Word, from_string
	from thompson import word
	from thompson.generators import Generators
	from thompson.orbits import dump_orbit_types
"""

__all__ = ["Automorphism"]

from collections import deque
from itertools import product
from io import StringIO

from .word import *
from .generators import Generators
from .orbits import *

#TODO. Check the assumption that the bases consist of simple words only (no lambdas).
#Expand until all the words are simple
class Automorphism:
	r"""Represents an automorphism of :math:`V_{n,r}` by specifying two bases. This class keeps track of the mapping between bases.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	:ivar domain: a :class:`generating set <thompson.generators.Generators>` of preimages 
	:ivar range: a :class:`generating set <thompson.generators.Generators>` of images.
	"""
	
	def __init__(self, arity, alphabet_size, domain, range):
		r"""Creates an automorphism, given the *arity* :math:`n` and *alphabet_size* :math:`r`. Two bases *domain* and *range* are given. The automorphism maps elements so that the given order of generators is preserved:
		
			.. math:: \text{domain}_i \mapsto \text{range}_i
		
		:raises ValueError: if the bases are of different sizes.
		:raises IndexError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis isn't actually a basis, i.e. is not a :meth:`free generating set <thompson.generators.Generators.is_free>` or does not :meth:`contract to the standard basis <thompson.generators.Generators.test_generates_algebra>`.
		"""
		
		#The boring checks
		if len(domain) != len(range):
			raise ValueError("Domain basis has {} elements, but range basis has {} elements.".format(
			  len(domain), len(range)))
		
		if not(arity == domain.arity == range.arity):
			raise IndexError("Arities do not match. Expected: {}, Domain: {}, Range: {}".format(
			  arity, domain.arity, range.arity))
		
		if not(alphabet_size == domain.alphabet_size == range.alphabet_size):
			raise IndexError("Alphabet sizes do not match. Expected: {}, Domain: {}, Range: {}".format(
			  alphabet_size, domain.alphabet_size, range.alphabet_size))
		
		#Check to see that the two generating sets given are free
		i, j = domain.test_free()
		if not(i == j == -1):
			raise ValueError("Domain is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		i, j = range.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		#Check to see that the two generating sets generate all of V_{n,r}
		missing = domain.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Domain does not generate V_{},{}. Missing elements are {}.".format(
			  arity, alphabet_size, [format(x) for x in missing]))
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{},{}. Missing elements are {}.".format(
			  arity, alphabet_size, [format(x) for x in missing]))
		
		#Before saving the domain and range, reduce them to remove any redundancy. This is like reducing tree pairs.
		#Is this really 100% necessary? Seem to remember the big lemma/algorithm for QNF needed it.
		#Yes, it is: how do you know that you've found the smallest possible nice basis if you haven't kept everything as small as possible throughout?
		Automorphism._reduce(domain, range)
		
		self.arity = arity
		self.alphabet_size = alphabet_size
		self.domain = domain
		self.range = range
		
		self._map = {}
		self._inv = {}
		for d, r in zip(self.domain, self.range):
			self._set_image(d, r)
		#Compute and cache the images of any element in X<A> above self.domain.
		for root in Generators.standard_basis(self.arity, self.alphabet_size):
			self._image_simple_above_domain(root)
			self._image_simple_above_domain(root, inverse=True)
		
		self._qnf_basis = None
		self._qnf_orbit_types = {}
		
		#Cache attributes
		#_map
		#_inv
		#_qnf_basis
		#_qnf_orbit_types
	
	@staticmethod
	def _reduce(domain, range):
		"""Contracts the domain generators whenever the corresponding contraction in range is possible. (This corresponds to reducing a tree pair diagram.)
			
			>>> from thompson.examples import cyclic_order_six as cyclic
			>>> #This is given by 6 generators, but after reduction consists of 5:
			>>> for d, r in zip(cyclic.domain, cyclic.range):
			... 	print(d, '->', r)
			x1 a1 a1 -> x1 a1 a1
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2
			x1 a2 a1 -> x1 a1 a2 a2 a1
			x1 a2 a2 -> x1 a1 a2 a1
			>>> #Swaps x1 and x2
			>>> domain = Generators(2, 2, ["x1", "x2 a1", "x2 a2"])
			>>> range  = Generators(2, 2, ["x2", "x1 a1", "x1 a2"])
			>>> Automorphism._reduce(domain, range)
			>>> for d, r in zip(domain, range):
			... 	print(d, '->', r)
			x1 -> x2
			x2 -> x1
		"""
		#similar to word._reduce and Generator.test_generates_algebra
		i = 0
		arity = domain.arity
		while i <= len(domain) - arity:
			d_pref = are_contractible(domain[i : i + arity])
			r_pref = are_contractible(range[i : i + arity])
			if d_pref and r_pref: #are both non-empty tuples
				domain[i : i + arity] = [Word(d_pref, arity, domain.alphabet_size)]
				range[ i : i + arity] = [Word(r_pref, arity, range.alphabet_size )]
				i -= (arity - 1) 
				i = max(i, 0)
			else:
				i += 1
	
	def _set_image(self, d, r, inverse=False):
		"""Stores the rule that phi(d) = r in the mapping dictionaries and ensure that d and r are both Words. Use inverse to specify that phi(r) = d."""
		assert isinstance(d, Word), repr(d)
		assert isinstance(r, Word), repr(r)
		if inverse:
			self._map[r] = d
			self._inv[d] = r
		else:
			self._map[d] = r
			self._inv[r] = d
	
	#Finding images of words
	#todo tests for inverse
	def image(self, key, inverse=False):
		"""Computes the image of a *key* under the given automorphism. The result is cached for further usage later. If *inverse* is True, the inverse of the current automorphism is used to map *key* instead.
		
		The *key* must be given as one of:
		
		- a sequence of integers (see the :mod:`~thompson.word` module), or
		- a :class:`~thompson.word.Word` instance.
		
		The input need not be in standard form. This method
		
		1. Converts *key* to standard form if necessary.
		2. Checks if the image of the standard form of *key* has been cached, and returns the image if so.
		3. Computes the image of *key* under the automorphism, then caches and returns the result.
		
			>>> from thompson.examples import example_4_25
			>>> #An element of the domain---just a lookup
			>>> print(example_4_25.image('x1 a1'))
			x1 a1 a1 a1
			>>> #A word below a the domain words
			>>> print(example_4_25.image('x1 a1 a2 a2'))
			x1 a1 a1 a1 a2 a2
			>>> #Above domain words---have to expand.
			>>> print(example_4_25.image('x1'))
			x1 a1 a1 a1 x1 a1 a1 a2 x1 a2 a2 x1 a1 a2 L x1 a2 a1 L L L
			>>> #Let's try some words not in standard form
			>>> print(example_4_25.image('x1 a1 a1 x1 a1 a2 L'))
			x1 a1 a1 a1
			>>> print(example_4_25.image('x1 a1 a1 x1 a1 a2 L a2 a1'))
			x1 a1 a1 a1 a2 a1
			>>> print(example_4_25.image('x1 a1 a1 x1 a2 a2 L'))
			x1 a1 a1 a1 a1 x1 a2 a2 x1 a1 a2 L x1 a2 a1 L L
		
		:rtype: a Word instance (which are always in standard form).
		"""
		dict = self._inv if inverse else self._map
		try:
			return dict[key]
		except KeyError:
			if isinstance(key, Word):
				word = key
			elif isinstance(key, (str, tuple)):
				word = Word(key, self.arity, self.alphabet_size)
		
		try:
			return dict[word]
		except KeyError:
			pass
		
		#1. First deal with the easy words (no lambdas).
		if word.is_simple():
			#During initialisation we compute the image of everything above the domain.
			#Thus *word* must be below the domain, so we remove alphas until we find something we know the image of.
			return self._image_simple_below_domain(word, inverse)
		#2. Words with lambdas in them are going to be a bit trickier.
		return self._image_of_lambda(word, inverse)
	
	def _image_simple_above_domain(self, word, inverse=False):
		r"""Finds the image of a *word* in standard form above *self.domain* by expansion.
		
		.. math:: w\phi = w\alpha\lambda\phi &= w\alpha_1     \dots w\alpha_n \lambda \phi \\
		                                     &= w\alpha_1\phi \dots w\alpha_n \phi \lambda
		
		Images are cached once computed.
		"""
		dict = self._inv if inverse else self._map
		try: 
			return dict[word]
		except KeyError:
			img_letters = _concat(self._image_simple_above_domain(child, inverse) for child in word.expand())
			#NOT in standard form.
			img_letters = standardise(img_letters, self.arity)
			image = Word(img_letters, self.arity, self.alphabet_size, preprocess=False)
			self._set_image(word, image, inverse)
			return image
	
	def _image_simple_below_domain(self, word, inverse=False):
		r"""This method takes a :meth:`simple <Word.is_simple>` *word* of the form :math:`d \alpha_i_1 \dots \alpha_i_m` where :math:`d` is the largest such word whose image has already been computed. The images of the words
		
		.. math:: d \alpha_i_1, d \alpha_i_1 \alpha_i_2, \dotsc, d\alpha_i_1 \alpha i_2 \dots \alpha i_m
		
		are then computed and cached. The final image in this list (i.e. that of the original *word* is returned).
		"""
		dict = self._inv if inverse else self._map
		i = 1
		while True:
			head, tail = word.rsplit(i)
			if head in dict:
				break
			i += 1
		head = Word(head, self.arity, self.alphabet_size, preprocess=False)
		image = dict[head] #is a word
		for _ in range(i):
			alpha, tail = tail[0], tail[1:]
			head = head.alpha(-alpha) #.alpha() returns a word
			image = image.alpha(-alpha) #.alpha() returns a word
			self._set_image(head, image, inverse)
		assert len(tail) == 0
		return image
		
	def _image_of_lambda(self, word, inverse=True):
		r"""This method accepts a *word* which ends in a lambda and extracts the arguments of the lambda. Next, the method calls :meth:`_get_image` on each argument, continuing the process recursively. Once we have computed/retreived all of the images we need, the images are concatenated and standardised, forming the image of *word*. This image is cached and returned.
		
		:raises ValueError: if the last letter in *word* is not a lambda.
		"""
		dict = self._inv if inverse else self._map
		subwords = lambda_arguments(word)
		letters = _concat(self.image(word, inverse) for word in subwords)
		letters = standardise(letters, self.arity)
		image = Word(letters, self.arity, self.alphabet_size, preprocess=False) #can skip validation
		self._set_image(word, image, inverse)
		return image
	
	def repeated_image(self, key, power):
		r"""If :math:`\psi` is the current automorphism, returns :math`\text{key}\psi^\text{power}`.
		
		:rtype: a :class:`~thompson.word.Word` instance.
		"""
		#TODO doctests.
		if power == 0:
			if not isinstance(key, Word):
				key = Word(key, self.arity, self.alphabet_size)
			return key
		inverse = power < 0
		power = abs(power)
		image = key
		for _ in range(power):
			image = self.image(image, inverse=inverse)
		return image
		
	#Printing
	def __str__(self):
		"""Printing an automorphism gives its arity, alphabet_size, and lists the images of its domain elements.
		
			>>> from thompson.examples import cyclic_order_six
			>>> print(cyclic_order_six)
			Automorphism of V_2,1 specified by 5 generators (after reduction):
			x1 a1 a1    -> x1 a1 a1
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2
			x1 a2 a1    -> x1 a1 a2 a2 a1
			x1 a2 a2    -> x1 a1 a2 a1
		"""
		output = StringIO()
		output.write("Automorphism of V_{},{} specified by {} generators (after reduction):".format(
		  self.arity, self.alphabet_size, len(self.domain)))
		max_len = 0
		for key in self.domain:
			max_len = max(max_len, len(str(key)))
		fmt = "\n{!s: <" + str(max_len) + "} -> {!s}"
		for key in sorted(self.domain):
			value = self.image(key)
			output.write(fmt.format(key, value))
		return output.getvalue()
	
	def dump_mapping(self, inverse=False, **kwargs):
		r"""A tool for debugging. Prints all the mapping rules :math:`w \mapsto \phi(w)` that have been stored. If *inverse* is True, the rules for the inverse are printed. Any keyword arguments *kwargs* are passed to the :func:`print <py3:print>` function."""
		dict = self._inv if inverse else self._map
		print('Automorphism of V_{},{} specified by {} generators (after reduction).'.format(
		  self.arity, self.alphabet_size, len(self.domain)), **kwargs)
		if inverse:
			print('Printing the inverse mapping.', **kwargs)
		print('The generators are marked by an asterisk.', **kwargs)
		max_len = 0
		for key in dict:
			max_len = max(max_len, len(str(key)))
		fmt = "{!s: <" + str(max_len) + "} -> {!s}"
		
		generators = self.range if inverse else self.domain
		
		for key in sorted(dict):
			prefix = " * " if key in generators else "   "
			value = self.image(key, inverse)
			print(prefix + fmt.format(key, value), **kwargs)
	
	#Operations on automorphisms
	def quasinormal_basis(self):
		r"""An implementation of Lemma 4.24.1. In [Hig]_ (section 9) Higman defines when an automorphism :math:`\phi` is in *quasinormal form* with respect to a given basis :math:`X`.  We return the basis :math:`X` w.r.t which the current automorphism is in quasinormal form.
		
		>>> from thompson.examples import *
		>>> example_4_5.quasinormal_basis()
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
		"""
		if self._qnf_basis is not None:
			return self._qnf_basis
		
		#TODO test me to high heaven.
		basis = self._minimal_expansion()
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
	
	def _minimal_expansion(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \cup W`. See example 4.25.
		
		>>> from thompson.examples import *
		>>> example_4_5._minimal_expansion()
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
		>>> example_4_11._minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2'])
		>>> example_4_12._minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2'])
		>>> example_4_25._minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
		>>> cyclic_order_six._minimal_expansion()
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2 a1', 'x1 a1 a2 a2', 'x1 a2'])
		"""
		#TODO arity 3 example, alphabet size > 1 example.
		basis = Generators.standard_basis(self.arity, self.alphabet_size)
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
		
		>>> from thompson.examples import *
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
		>>> basis = example_4_12._minimal_expansion()
		>>> dump_orbit_types(example_4_12, basis)
		x1 a1: Incomplete orbit
		x1 a2: Incomplete orbit
		with respect to the basis [x1 a1, x1 a2]
		>>> basis.expand(0)
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2', 'x1 a2'])
		>>> dump_orbit_types(example_4_12, basis)
		x1 a1 a1: Bi-infinite orbit containing [x1 a2] a2
		x1 a1 a2: Bi-infinite orbit containing [x1 a2] a1
		x1 a2: Incomplete orbit
		with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2]
		>>> basis.expand(2)
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
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
	
	def share_orbit(self, u, v):
		r"""Uses lemma 4.25 part (2) to determine if :math:`u` and :math:`v` are in the same orbit of the current automorphism :math:`\psi`. Specifically, does there exist an integer :math:`n` such that :math:`u\psi^n = v`?
		
			>>> u  = Word('x a2 a2 a1 a1 a2', 2, 1)
			>>> v1 = Word('x a1 a2', 2, 1)
			>>> v2 = Word('x a1 a1 a2', 2, 1)
			>>> print(example_4_25.share_orbit(u, v1))
			None
			>>> example_4_25.share_orbit(u, v2)
			3
			
		
		:returns: An integer :math:`n` if it exists; otherwise ``None``. In particular if :math:`u = v` then this method returns ``0``. 
		"""
		# print('u = {}\nv = {}'.format(u, v))
		if u == v:
			return 0
		basis = self.quasinormal_basis()
		# print('QNF basis is', basis)
		if not basis.is_above(u):
			# u = u._as_contraction()
			raise ValueError('todo')
		if not basis.is_above(v):
			raise ValueError('todo')
			# v = v._as_contraction()
			
		if not (u.is_simple() and v.is_simple()):
			raise ValueError('todo')
			'''depth = max(u.max_depth_to(basis), v.max_depth_to(basis))
			alphas = range(-1, -self.arity, -1)
			solution_set = TODO_ALL_INTEGERS #TODO
			# For all strings of length *depth* \Gamma made only from alphas:
			for tail in product(alphas, repeat=depth):
				u_desc = u.extend(tail)
				v_desc = v.extend(tail)
				solution &= share_orbit(u_desc, v_desc)
				if solution_set is not EMPTY: #TODO
					return solution_set
			return solution_set'''
		
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
				return None
			
			period = u_head_type.data
			image = u
			# print('Testing the orbit of u.')
			for i in range(1, period):
				image = self.image(image)
				if image == v:
					# print('Found a match.')
					return i
			# print('No match found.')
			return None
		
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
				return u_shift + i - v_shift
		return None
		
	def _preprocess(self, head, tail):
		r"""Takes a pair :math:`(y, \Gamma)` below the quasi-normal basis :math:`X` and returns a triple :math:`(\widetilde{y}, \widetilde{\Gamma}, k) where
		* The orbit of :math:`\widetilde{y}` is semi-infinite;
		* :math:`\widetilde{y}\widetilde{\Gamma}` is in the same orbit as :math:`y\Gamma`.
		* :math:`\widetilde{\Gamma}` does not start with the characteristic multiplier of :math:`\widetilde{y}`.
		
			>>> example_4_25.quasinormal_basis()
			Generators(2, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
			>>> example_4_25._preprocess(Word('x a2 a2', 2, 1), from_string('a1 a1 a2'))
			(1, Word('x1 a2 a2', 2, 1), (-2,))
		"""
		head_type = self._qnf_orbit_types[head]
		if head_type.is_type('C'):
			print('infinite orbit: have to find type B child')
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
		"""Makes a copy of *tail* and removes as many copies of *multiplier* from the start of the copy as possible. Returns the copy.
			
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

#TODO. Compose and invert automorphisms. Basically move all the functionality from TreePair over to this class and ditch trree pair.
#TODO method to decide if the automorphism is in (the equivalent of) F, T, or V.
#TODO the named elements A, B, C, X_n of Thompson's V.

