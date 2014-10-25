r"""The algebra of :mod:`words <thompson.word>` has its own structure, and just like groups, rings etc. there exist maps (homomorphisms) which preserve this structure. In our specific case, a homomorphism :math:`\phi: V_{n,r} \to V_{n,s}` is a function satisfying

.. math:: w \alpha_i \phi = w \phi \alpha_i \qquad \text{and} \qquad
          w_1 \dots w_n \lambda \phi = w_1 \phi \dots w_n \phi \lambda.

In other words, a homomorphism is a function which 'commutes' with the algebra operations :math:`\alpha_i` and :math:`\lambda`.

As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. 

This module works with automorphisms of :math:`V_{n,r}(\boldsymbol{x})`. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.
 
.. testsetup::
	
	from thompson.word import Word, from_string
	from thompson import word
	from thompson.generators import Generators
	from thompson.automorphism import *
	from thompson.orbits import dump_orbit_types
	from thompson.examples import *
"""

__all__ = ["Isomorphism", "Automorphism"]

from collections import defaultdict, deque
from io import StringIO
from itertools import product

from .word import *
from .generators import Generators
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

class Isomorphism:
	r"""Represents an ismorphism of algebras :math:`\psi : V_{n,r} \to V_{n,s}` as a bijection between bases.
	:ivar domain_words: a :class:`basis <thompson.generators.Generators>` of preimages 
	:ivar range_words:  a :class:`basis <thompson.generators.Generators>` of images.
	
	:raises ValueError: if the bases are of different sizes.
	:raises TypeError: if the algebras described by *domain* and *range* are not isomorphic.
	"""
	def __init__(self, domain, range):
		#The boring checks
		if len(domain) != len(range):
			raise ValueError("Domain basis has {} elements, but range basis has {} elements.".format(
			  len(domain), len(range)))
		
		if not domain.signature.is_isomorphic_to(range.signature):
			raise TypeError("Domain signature {} is not isomorphic to range signature {}.".format(
			  domain.signature, range.signature))
		
		#Expand any non-simple words
		Isomorphism._expand(domain, range)
		#Before saving the domain and range, reduce them to remove any redundancy. This is like reducing tree pairs.
		#How do you know that you've found the smallest possible nice basis if you haven't kept everything as small as possible throughout?
		Isomorphism._reduce(domain, range)
		
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
			raise ValueError("Domain does not generate V_{}. Missing elements are {}.".format(
			  domain.signature, [format(x) for x in missing]))
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{}. Missing elements are {}.".format(
			  range.signature, [format(x) for x in missing]))
		
		self.domain = domain
		self.range = range
		
		#Mapping caches
		self._map = {}
		self._inv = {}
		for d, r in zip(self.domain, self.range):
			self._set_image(d, r)
			
		#Compute and cache the images of any simple word above self.domain.
		#Setup the inverse map
		for root in Generators.standard_basis(domain.signature):
			self._image_simple_above_domain(root)
			self._image_simple_above_domain(root, inverse=True)
	
	@staticmethod
	def _expand(domain, range):
		"""Expands the pair of generating sets where neccessary until all words in both sets are simple."""
		#todo doctest and string
		i = 0
		while i < len(domain):
			if not (domain[i].is_simple() and range[i].is_simple()):
				domain.expand(i)
				range.expand(i)
			else:
				i += 1
	
	@staticmethod
	def _reduce(domain, range):
		"""Contracts the domain generators whenever the corresponding contraction in range is possible. This corresponds to reducing a tree pair diagram.
			
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
			>>> domain = Generators((2, 2), ["x1", "x2 a1", "x2 a2"])
			>>> range  = Generators((2, 2), ["x2", "x1 a1", "x1 a2"])
			>>> Automorphism._reduce(domain, range)
			>>> for d, r in zip(domain, range):
			... 	print(d, '->', r)
			x1 -> x2
			x2 -> x1
		"""
		#similar to word._reduce and Generator.test_generates_algebra
		#This method ensures that self.domain.minimal_expansion(self) == self.domain after initialisation.
		i = 0
		arity = domain.signature.arity
		while i <= len(domain) - arity:
			d_pref = are_contractible(domain[i : i + arity])
			r_pref = are_contractible(range[ i : i + arity])
			if d_pref and r_pref: #are both non-empty tuples
				domain[i : i + arity] = [Word(d_pref, domain.signature)]
				range[ i : i + arity] = [Word(r_pref, range.signature)]
				i -= (arity - 1) 
				i = max(i, 0)
			else:
				i += 1
	
	def _set_image(self, d, r, inverse=False):
		"""Stores the rule that phi(d) = r in the mapping dictionaries and ensure that d and r are both Words. Use inverse to specify that phi(r) = d."""
		if inverse:
			d, r = r, d
		assert d in self.domain.signature, repr(d)
		assert r in self.range.signature, repr(r)
		self._map[d] = r
		self._inv[r] = d
	
	#Simple operations on isomorphisms
	def __eq__(self, other):
		return self.domain == other.domain and self.range == other.range
	
	def __mul__(self, other): #self * other is used for the (backwards) composition self then other
		#TODO doctest and string
		if not isinstance(other, Isomorphism):
			return NotImplemented
		if self.range.signature != other.domain.signature:
			raise TypeError("Signatures do not match.")
		
		range = Generators(other.range.signature)
		for d, word in zip(self.domain, self.range):
			range.append(other.image(word))
		return Isomorphism(self.domain, range)
		#TODO will need to be think about when to return automorphisms here
	
	#Finding images of words
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
		
		Examples of finding inverse images:
		
			>>> print(example_4_25.image('x1 a2 a2', inverse=True))
			x1 a2 a2 a1 a1
			>>> print(example_4_25.image('x1 a1 a1 a2 a2 a1', inverse=True))
			x1 a2 a1 a2 a1
			>>> print(example_4_25.image('x a2', inverse=True))
			x1 a2 a2 a2 x1 a2 a2 a1 a1 L
			>>> print(example_4_25.image('x a2 a2 x a1 a2 L', inverse=True))
			x1 a2 a2 a1
		
		:rtype: a :class:`~thompson.word.Word` instance (which are always in standard form).
		"""
		dict = self._inv if inverse else self._map
		signature_input = self.range.signature if inverse else self.domain.signature
		try:
			return dict[key]
		except KeyError:
			if isinstance(key, Word):
				if key.signature != signature_input:
					raise TypeError("Signature {} of {} does not match {}".format(
					  key.signature, key, signature_input))
				word = key
			elif isinstance(key, (str, tuple)):
				word = Word(key, signature_input)
		
		try:
			return dict[word]
		except KeyError:
			pass
		
		#1. First deal with the easy words (no lambdas).
		#TODO. I have a feeling that this ought to be ``if self.domain.is_above(word) or something similar...
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
		signature_output = self.domain.signature if inverse else self.range.signature
		try: 
			return dict[word]
		except KeyError:
			img_letters = _concat(self._image_simple_above_domain(child, inverse) for child in word.expand())
			#NOT in standard form.
			img_letters = standardise(img_letters, word.signature)
			image = Word(img_letters, signature_output, preprocess=False)
			self._set_image(word, image, inverse)
			return image
	
	def _image_simple_below_domain(self, word, inverse=False):
		r"""This method takes a :meth:`simple <Word.is_simple>` *word* of the form :math:`d \alpha_i_1 \dots \alpha_i_m` where :math:`d` is the largest such word whose image has already been computed. The images of the words
		
		.. math:: d \alpha_i_1, d \alpha_i_1 \alpha_i_2, \dotsc, d\alpha_i_1 \alpha i_2 \dots \alpha i_m
		
		are then computed and cached. The final image in this list (i.e. that of the original *word* is returned).
		"""
		dict = self._inv if inverse else self._map
		signature_input = self.range.signature if inverse else self.domain.signature
		signature_output = self.domain.signature if inverse else self.range.signature
		i = 1
		while True:
			head, tail = word.rsplit(i)
			if head in dict:
				break
			i += 1
		head = Word(head, signature_input, preprocess=False)
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
		signature_output = self.domain.signature if inverse else self.range.signature
		subwords = lambda_arguments(word)
		letters = _concat(self.image(word, inverse) for word in subwords)
		letters = standardise(letters, self.signature)
		image = Word(letters, signature_output, preprocess=False) #can skip validation
		self._set_image(word, image, inverse)
		return image
	
	#Printing
	def _string_header(self):
		return "{}: V{} -> V{} specified by {} generators (after reduction).".format(
		  type(self).__name__, self.domain.signature, self.range.signature, len(self.domain))
	
	def __str__(self):
		"""Printing an automorphism gives its arity, alphabet_size, and lists the images of its domain elements.
		
			>>> from thompson.examples import cyclic_order_six
			>>> print(cyclic_order_six)
			Automorphism: V(2, 1) -> V(2, 1) specified by 5 generators (after reduction).
			x1 a1 a1    -> x1 a1 a1
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2
			x1 a2 a1    -> x1 a1 a2 a2 a1
			x1 a2 a2    -> x1 a1 a2 a1
		"""
		output = StringIO()
		output.write(self._string_header())
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
		print(self._string_header(), **kwargs)
		if inverse:
			print('Printing the **inverse** mapping.', **kwargs)
		
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

class Automorphism(Isomorphism):
	r"""An automorphism is an isomorphism from :math:`V_{n,r}` to itself. The automorphism group :math:`\mathrm{Aut}(V_{n,r})` is labelled :math:`G_{n,r}`.
	"""
	def __init__(self, domain, range):
		r"""Creates an automorphism, given the *arity* :math:`n` and *alphabet_size* :math:`r`. Two bases *domain* and *range* are also given. The automorphism maps elements so that the given order of generators is preserved:
		
			.. math:: \text{domain}_{\,i} \mapsto \text{range}_{\,i}
		
		:raises IndexError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis :meth:`isn't actually a basis <thompson.generators.Generators.is_basis>`.
		"""
		if domain.signature != range.signature:
			raise TypeError('Domain {} and range {} have different signatures.'.format(
			  domain.signature, range.signature))
		
		super().__init__(domain, range)
		self.signature = domain.signature
		
		#Cache attributes
		self._qnf_basis = None
		self._qnf_orbit_types = {}
	
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
			image = self.image(image, inverse=inverse)
		return image
	
	#Generating the quasinormal basis
	def quasinormal_basis(self):
		r"""We say that :math:`\phi` is *in semi-normal form* with respect to the basis :math:`X` if no element of :math:`X` lies in an incomplete :math:`X`-component of a :math:`\phi` orbit. See the :mod:`~thompson.orbits` module for more details.
		
		There is a minimal such basis, :math:`X_\phi` say, and we say that :math:`\phi` is *in quasi-normal form* with respect to :math:`X_\phi`. This method determines and returns the basis :math:`X_\phi` where :math:`\phi` denotes the current automorphism. The result is cached so that further calls to this method perform no additional computation. 
		
			>>> example_4_5.quasinormal_basis()
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> alphabet_size_two.quasinormal_basis()
			Generators((3, 2), ['x1 a1', 'x1 a2', 'x1 a3', 'x2'])
		
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
		
		.. seealso:: This is an implementation of Algorithm 5.6 in the paper. It depends on Algorithms 5.13 and 5.27. See also :meth:`PeriodicFactor.test_conjugate_to` and :meth:`InfiniteFactor.test_conjugate_to`.
		"""
		#todo doctests
		#TODO broken links.
		#0. Check that both automorphisms belong to the same G_{n,r}.
		if not(self.arity == other.arity
		  and self.alphabet_size == other.alphabet_size):
			raise ValueError('Both automorphisms must have the same arity and alphabet size.')
		
		#1. Before going any further, check that the number of periodic and infinite elements are compatible.
		result = self._check_parition_sizes(other)
		if result is None:
			return None
		pure_periodic, pure_infinite, s_qnf_p, s_qnf_i, o_qnf_p, o_qnf_i = result
		
		#4. If necessary, test the periodic factors.
		if not pure_infinite:
			#2, 3. Extract the periodic factor.
			s_p = self.free_factor(s_qnf_p, infinite=False)
			o_p = self.free_factor(o_qnf_p, infinite=False)
			
			rho_p = s_p.test_conjugate_to(o_p)
			if rho_p is None:
				return None
			#TODO rewrite rho_p in terms of the original algebra
			if pure_periodic:
				return rho_p
		
		#Step 5. If necessary, test the infinite factors.
		if not pure_periodic:
			#2, 3. Extract the infinite factor.
			s_i = self.free_factor(s_qnf_i, infinite=True)
			o_i = self.free_factor(o_qnf_i, infinite=True)
			
			rho_i = s_i.test_conjugate_to(o_i)
			if rho_i is None:
				return None
			#TODO rewrite rho_i in terms of the original algebra
			if pure_infinite:
				return rho_i
		
		#Step 6. If we've got this far, we have a periodic/infinite mixture.
		#Combine the rewritten conjugators into one conjugator in G_n,r
		return rho_p * rho_i
	
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
		modulus = self.arity - 1
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
			Using {y1} to denote the root letters of this automorphism.
			PeriodicFactor: V(2, 1) -> V(2, 1) specified by 2 generators (after reduction).
			y1 a1 -> y1 a2
			y1 a2 -> y1 a1
			This is embedded in a parent automorphism by the following rules.
			y1 a1 ~~> x1 a1 a2 a1
			y1 a2 ~~> x1 a1 a2 a2
		
		.. doctest::
			:hide:
			
			>>> # alphabet_size_two example
			>>> qnb = alphabet_size_two.quasinormal_basis()
			>>> p, i = alphabet_size_two._partition_basis(qnb)
			>>> print(alphabet_size_two.free_factor(p, infinite=False))
			Using {y1} to denote the root letters of this automorphism.
			PeriodicFactor: V(3, 1) -> V(3, 1) specified by 1 generators (after reduction).
			y1 -> y1
			This is embedded in a parent automorphism by the following rules.
			y1 ~~> x1 a1
			>>> print(alphabet_size_two.free_factor(i, infinite=True))
			Using {y1} to denote the root letters of this automorphism.
			InfiniteFactor: V(3, 1) -> V(3, 1) specified by 5 generators (after reduction).
			y1 a1    -> y1 a1 a3
			y1 a2    -> y1 a3
			y1 a3 a1 -> y1 a2
			y1 a3 a2 -> y1 a1 a2
			y1 a3 a3 -> y1 a1 a1
			This is embedded in a parent automorphism by the following rules.
			y1 a1    ~~> x1 a2
			y1 a2    ~~> x1 a3
			y1 a3 a1 ~~> x2 a1
			y1 a3 a2 ~~> x2 a2
			y1 a3 a3 ~~> x2 a3
		"""
		#TODO more doctests
		if len(generators) == 0:
			raise ValueError('Must provide at least one generator.')
		expansion = generators.minimal_expansion_for(self)
		modulus = self.signature.arity - 1
		alphabet_size = modulo_non_zero(len(generators), modulus)
		
		images = Generators.standard_basis((self.signature.arity, alphabet_size))
		images.expand_to_size(len(generators))
		domain, range = self._rewrite_mapping(expansion, generators, images)
		
		type = InfiniteFactor if infinite else PeriodicFactor
		return type(domain, range, expansion)
	
	def _rewrite_mapping(self, words, basis, img_basis):
		r"""Given a *basis*, a set of *words* below that *basis*, and a bijective image *img_basis* of *basis*, this method rewrites the list of rules ``w -> self[w] for w in words`` in terms of the new *img_basis*.
		
		:returns: a pair of lists *domain, range* where *domain* and *range* are both below *img_basis*.
		"""
		domain = Generators(img_basis.signature)
		range = Generators(img_basis.signature)
		
		for word in words:
			preimage = self._rewrite(word, basis, img_basis)
			domain.append(preimage)
			image = self._rewrite(self.image(word), basis, img_basis)
			range.append(image)
		
		return domain, range
	
	def _rewrite(self, word, basis, img_basis):
		r"""Suppose that we have a *word* which is below a given *basis*. Suppose we have a bijection between *basis* and some image *img_basis*. Then we can rewrite *word* as as descendant of *img_basis* in a manner which is compatible with this bijection.
		"""
		result = basis.test_above(word, return_index=True)
		if result is None:
			raise ValueError('The word {} is not below the basis {}'.format(
			  word, basis))
		index, tail = result
		return img_basis[index].extend(tail)

#TODO? Compose and invert automorphisms. Basically move all the functionality from TreePair over to this class and ditch trree pair.
#TODO method to decide if the automorphism is in (the equivalent of) F, T, or V.
#TODO the named elements A, B, C, X_n of Thompson's V.

class AutomorphismFactor(Automorphism):
	"""An automorphism derived from a larger parent automorphism.
	
	:ivar domain_preimage: a list of words which tells us where each domain word originally came from.
	"""
	def __init__(self, domain, range, domain_preimage):
		"""In addition to creating an automorphism, we store the *domain_preimage* to keep track of the mapping between the words in this automorphism's algebra and the words manipulated by the parent automorphism.
		
		:raises ValueError: if ``len(domain_preimage) != len(domain)``.
		
		.. seealso:: :meth:`Automorphism.__init__`
		"""
		#todo example
		if not len(domain_preimage) == len(domain):
			raise ValueError('Each domain word must be associated to a word from the parent algebra. Found {} associations but domain size is {}'.format(
			  len(domain_preimage), len(domain)))
		super().__init__(domain, range)
		self.domain_preimage = domain_preimage
	
	def __str__(self):
		output = StringIO()
		if self.signature.alphabet_size == 1:
			roots = '{y1}'
		else:
			roots = '{{y1, ..., y{}}}'.format(self.alphabet_size)
		output.write("Using {} to denote the root letters of this automorphism.\n".format(
		  roots))
		output.write(super().__str__().replace('x', 'y'))
		output.write("\nThis is embedded in a parent automorphism by the following rules.")
		
		max_len = 0
		for key in self.domain:
			max_len = max(max_len, len(str(key)))
		fmt = "\n{: <" + str(max_len) + "} ~~> {!s}"
		
		for domain_word, preimage in zip(self.domain, self.domain_preimage):
			output.write(fmt.format(
			  str(domain_word).replace('x', 'y'), preimage))
		return output.getvalue()

class PeriodicFactor(AutomorphismFactor):
	r"""A purely periodic free factor which has been extracted from another component.
	
		>>> print(example_5_9_p)
		Using {y1} to denote the root letters of this automorphism.
		PeriodicFactor: V(2, 1) -> V(2, 1) specified by 7 generators (after reduction).
		y1 a1 a1    -> y1 a1 a2 a1
		y1 a1 a2 a1 -> y1 a1 a2 a2
		y1 a1 a2 a2 -> y1 a1 a1
		y1 a2 a1 a1 -> y1 a2 a1 a2
		y1 a2 a1 a2 -> y1 a2 a1 a1
		y1 a2 a2 a1 -> y1 a2 a2 a2
		y1 a2 a2 a2 -> y1 a2 a2 a1
		This is embedded in a parent automorphism by the following rules.
		y1 a1 a1    ~~> x1 a1 a1 a1
		y1 a1 a2 a1 ~~> x1 a1 a1 a2
		y1 a1 a2 a2 ~~> x1 a1 a2
		y1 a2 a1 a1 ~~> x1 a2 a1 a1
		y1 a2 a1 a2 ~~> x1 a2 a1 a2
		y1 a2 a2 a1 ~~> x1 a2 a2 a1
		y1 a2 a2 a2 ~~> x1 a2 a2 a2
		>>> sorted(example_5_9_p.cycle_type)
		[2, 3]
		>>> from pprint import pprint
		>>> #Two orbits of size 2, one orbit of size 3
		>>> pprint(example_5_9_p.multiplicity)
		{2: 2, 3: 1}
		
	
	:ivar multiplicity: a mapping :math:`d \mapsto m_\phi(d, X_\phi)` where :math:`\phi` is the current automorphism and :math:`X_\phi` is the :meth:`quasi-normal basis <thompson.automorphism.Automorphism.quasinormal_basis>` for :math:`\phi`.
	:ivar cycle_type: the set :math:`\{d \in \mathbb{N} : \text{$\exists$ an orbit of length $d$.}\}`
	"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		#See discussion before lemma 5.3
		assert self.quasinormal_basis() == self.domain
		assert self.quasinormal_basis().minimal_expansion_for(self) == self.domain
		
		#Attributes:
		#.cycle_types
		#.multiplicities
		#Need to have constructed the QNF basis here.
		#TODO: can we just pass on this information (QNF basis, QNF orbit types) from the parent automorphism?
		self._setup_cycle_type()
		self._setup_multiplicities()
	
	def _setup_cycle_type(self):
		#see definition 5.8
		self.cycle_type = {type.data for type in self._qnf_orbit_types.values()}
	
	def _setup_multiplicities(self):
		counts = defaultdict(int)
		for orbit in self._qnf_orbit_types.values():
			orbit_size = orbit.data
			counts[orbit_size] += 1
		for orbit_size, occurances in counts.items():
			assert occurances % orbit_size == 0
			counts[orbit_size] = occurances // orbit_size
		self.multiplicity = dict(counts)
	
	def enumerate_orbits(self, basis):
		r"""Enumerates the periodic orbits of the current automorphism's quasinormal_basis. Returns a dictionary *orbits_by_size*. Each value ``orbits_by_size[d]`` is a list of the orbits of size *d*. Orbits themselves are represented as lists of :class:`Words <thompson.word.Word>`.
		
			>>> def display_orbits(orbits_by_size):
			... 	for key in sorted(orbits_by_size):
			... 		print('Orbits of length', key)
			... 		for list in orbits_by_size[key]:
			... 			print('...', *list, sep=' -> ', end=' -> ...\n')
			>>> orbits_by_size = example_5_9_p.enumerate_orbits(example_5_9_p.quasinormal_basis())
			>>> display_orbits(orbits_by_size)
			Orbits of length 2
			... -> x1 a2 a1 a1 -> x1 a2 a1 a2 -> ...
			... -> x1 a2 a2 a1 -> x1 a2 a2 a2 -> ...
			Orbits of length 3
			... -> x1 a1 a1 -> x1 a1 a2 a1 -> x1 a1 a2 a2 -> ...
		"""
		#TODO doctest
		#TODO check this is deterministic
		orbits_by_size = defaultdict(deque)
		already_seen = set()
		for gen in basis:
			if gen in already_seen:
				continue
			type, images = self._orbit_type(gen, basis)
			length = type.data
			images = [images[i] for i in range(length)]
			already_seen.update(images)
			orbits_by_size[length].append(images)
		return dict(orbits_by_size)
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely periodic automorphisms are periodic by examining their orbits based on their size.
		
			>>> psi_p = example_5_12_psi_p; phi_p = example_5_12_phi_p
			>>> rho_p = example_5_12_psi_p.test_conjugate_to(example_5_12_phi_p)
			>>> print(rho_p)
			Automorphism: V(2, 1) -> V(2, 1) specified by 6 generators (after reduction).
			x1 a1 a1    -> x1 a1 a2
			x1 a1 a2    -> x1 a2 a2
			x1 a2 a1 a1 -> x1 a1 a1 a1
			x1 a2 a1 a2 -> x1 a2 a1 a1
			x1 a2 a2 a1 -> x1 a1 a1 a2
			x1 a2 a2 a2 -> x1 a2 a1 a2
			>>> rho_p * phi_p == psi_p * rho_p
			True
		
		.. seealso:: This implements algorithm 5.13 of the paper---see section 5.3.
		"""
		# todo docstring and doctest
		if not isinstance(other, PeriodicFactor):
			raise TypeError('Other automorphism must be a PeriodicFactor.')
		
		# 1. The quasi-normal bases are constructed in initialisation.
		# Note that the QNF basis should be just the domain. TODO Checkme
		
		# 2. Check that the cycle types are the same.
		if self.cycle_type != other.cycle_type:
			return None
		
		#3. Check that the multiplicites are congruent.
		modulus = self.signature.arity - 1
		for d in self.cycle_type:
			if self.multiplicity[d] % modulus != other.multiplicity[d] % modulus:
				return None
		
		# 4. Expand bases until the orbits multiplicites are the same
		s_orbits_of_size = self.enumerate_orbits(self.quasinormal_basis())
		o_orbits_of_size = other.enumerate_orbits(other.quasinormal_basis())

		for d in self.cycle_type:
			assert len(s_orbits_of_size[d]) == self.multiplicity[d]
			assert len(o_orbits_of_size[d]) == other.multiplicity[d]
			
			expansions_needed = (self.multiplicity[d] - other.multiplicity[d]) // modulus
			
			if expansions_needed > 0:
				other.expand_orbits(o_orbits_of_size[d], expansions_needed)
			elif expansions_needed < 0:
				self.expand_orbits(s_orbits_of_size[d], -expansions_needed)
			
			assert len(s_orbits_of_size[d]) == len(o_orbits_of_size[d])
		
		domain = Generators(self.signature)
		range  = Generators(self.signature)
		for d in self.cycle_type:
			for s_orbit, o_orbit in zip(s_orbits_of_size[d], o_orbits_of_size[d]):
				for s_word, o_word in zip(s_orbit, o_orbit):
					domain.append(s_word)
					range.append(o_word)
		
		rho = Automorphism(domain, range)
		return rho
		#TODO translate back to original setting
	
	@staticmethod
	def expand_orbits(deque, num_expansions):
		r"""Takes a *deque* whose elements are a list of words. The following process is repeated *num_expansions* times.
		
		1. Pop an orbit :math:`\mathcal{O} = [o_1, \dotsc, o_d]` from the left of the deque.
		2. For :math:`i = 1, \dotsc, n` (where :math:`n` is the arity):
			a. Compute the new orbit :math:`\mathcal{O_i} =[o_1\alpha_i, \dotsc, o_d\alpha_i]`
			b. Append this orbit to the right of *deque*.
		
		Once complete, the number of orbits in *deque* has increased by :math:`\text{(arity -1) $\times$ num_expansions}`
		"""
		for _ in range(num_expansions):
			orbit = deque.popleft()
			arity = orbit[0].signature.arity
			for i in range(1, arity+1):
				new_orbit = [w.alpha(i) for w in orbit]
				deque.append(new_orbit)

class InfiniteFactor(AutomorphismFactor):
	def test_conjugate_to(self, other):
		raise NotImplementedError()