r"""The algebra of :mod:`words <thompson.word>` has its own structure, and just like groups, rings etc. there exist maps (homomorphisms) which preserve this structure. In our specific case, a homomorphism :math:`\phi: V_{n,r} \to V_{n,s}` is a function satisfying

.. math:: w \alpha_i \phi = w \phi \alpha_i \qquad \text{and} \qquad
          w_1 \dots w_n \lambda \phi = w_1 \phi \dots w_n \phi \lambda.

As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. 

This module works with automorphisms of :math:`V_{n,r}(\boldsymbol{x})`. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.
 
.. testsetup::
	
	from thompson.automorphism import *
	from thompson.word import Word
	from thompson.generators import Generators
	from thompson.orbits import dump_orbit_types
"""

__all__ = ["Automorphism"]#, "OrbitType", "orbit_types"]

from collections import deque
from itertools import chain, product
from io import StringIO

from .word import *
from .generators import Generators
from .orbits import *
from .full_tree import FullTree

def _concat(words):
	"""Takes an iterable *words* which yields lists of integers representing words. Returns a tuple containing all the *words* concatenated together, with a zero (lambda) added on the end."""
	return tuple(chain.from_iterable(words)) + (0,)

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
			type = self._orbit_type(basis[i], basis)
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
		"""Returns the orbit type of *y* with respect to the given *basis*.
		
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
		# print('Forward Orbit for', y)
		right_result = self._test_semi_infinite(y, basis, backward=False)
		# print('right_result:', right_result)
		# print('Backward orbit for', y)
		left_result  = self._test_semi_infinite(y, basis, backward=True)
		# print('left_result:', left_result)
		
		#A bunch of sanity checks. 
		#Suppose we fall out of the tree on the right.
		if right_result is False:
			if left_result is False:
				return OrbitType.incomplete
			
			assert left_result.type != OrbitType._complete_finite
			assert left_result.type != OrbitType._complete_infinite
			return left_result #Must have left semi-infinite remaining
		
		#Else right_result is an OrbitType instance
		if right_result.type == OrbitType._complete_infinite:
			assert left_result.type == OrbitType._complete_infinite
		
		if right_result.type == OrbitType._right_semi_infinite:
			assert left_result is False
		
		if right_result.type == OrbitType._complete_finite:
			assert left_result.type == OrbitType._complete_finite
			assert left_result.data == right_result.data #check they have the same period
		
		return right_result
		
	
	#TODO maybe this should return characteristics as well
	def _test_semi_infinite(self, y, basis, backward=False):
		"""Computes the orbit type of *y* with respect to *basis* in the forward direction. (Use ``backward=True`` to go backwards.)
		
		- If an incomplete orbit is found, returns False. Note we cannot conclude the orbit is incomplete without examining the other direction.
		- If a periodic orbit is found, returns an OrbitType containing the period.
		- If an infinite orbit is found:
			- if a characteristic was found, returns the appropriate semi-infinite OrbitType.
			- if no characteristic was found, returns OrbitType.complete_infinite.
		
		
		if a periodic orbit in X<A> is found; otherwise returns True or false to say if the orbit in X<A> is semi-infinite in the given direction."""
		image = y
		images = [y]
		m = 0
		
		while True:
			m += 1
			image = self.image(image, inverse=backward) #Compute y\phi^{\pm m}
			#1. Did we fall out of X<A>?
			if not basis.is_above(image):
				return False 
			
			#2. Is this a periodic orbit?
			if image == y:
				return OrbitType.periodic(m)
			
			#3. Is this an infinite orbit?
			#Find the generator *prefix* above the current image.
			prefix, image_tail = basis.test_above(image)
			for ell, previous in enumerate(images):
				previous_tail = prefix.test_above(previous)
				if previous_tail is not None:
					if previous is y:
						characteristic = self._extract_characteristic(prefix, m, ell, backward)
						return OrbitType.semi_infinite(characteristic, backward)
					else:
						type_b_data = prefix, previous_tail
						return OrbitType.complete_infinite(type_b_data)
						#At this point we know what the type of prefix is, too.
			#Otherwise, continue with the next power of \psi
			images.append(image)
	
	def _extract_characteristic(self, prefix, m, ell, backward):
		r"""Uses Lemma 4.14(B) to compute the characteristic :math:`(m, \Gamma)` of the orbit containing ``prefix + image_tail``.
		
		..warning:: This method is 
		"""
		#todo doctests
		i = m - ell
		if backward:
			i = -i
		image = self.repeated_image(prefix, i)
		# print(image)
		tail = prefix.test_above(image)
		assert tail is not None
		return i, tail
	
	def share_orbit(self, u, v):
		r"""Uses lemma 4.25 part (2) to determine if :math:`u` and :math`v` are in the same orbit of the current automorphism :math:\psi`.
		
		:returns: the set of powers :math:`m` such that :math:`u\psi^m = v`.
		"""
		if u == v:
			return TODO_ALL_INTEGERS  #TODO
		basis = self.quasinormal_basis()
		
		if not basis.is_above(u):
			u = u._as_contraction()
		if not basis.is_above(v):
			v = v._as_contraction()
			
		if not (u.is_simple() and v.is_simple()):
			depth = max(u.max_depth_to(basis), v.max_depth_to(basis))
			alphas = range(-1, -self.arity, -1)
			solution_set = TODO_ALL_INTEGERS #TODO
			#For all strings of length *depth* \Gamma made only from alphas:
			for tail in product(alphas, repeat=depth):
				u_desc = u.descend(tail)
				v_desc = v.descend(tail)
				solution_set &= share_orbit(u_desc, v_desc)
				if solution_set is EMPTY: #TODO
					return solution_set
			return solution_set
		
		#Now we're dealing with simple words below the basis.
		u_head, u_tail = basis.test_above(u) #y, gamma
		v_head, v_tail = basis.test_above(v) #z, delta
		
		
		
		#replace u v with semi-infinite equivalents such that u, v share orbit iff their equivalents do
		#remove copies of characteristic multipliers from u, v
		# proceed with the rest of the lemma?

#TODO. Compose and invert automorphisms. Basically move all the functionality from TreePair over to this class and ditch trree pair.
#TODO method to decide if the automorphism is in (the equivalent of) F, T, or V.
#TODO the named elements A, B, C, X_n of Thompson's V.

