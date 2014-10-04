r"""The algebra of :mod:`words <thompson.word>` has its own structure, and just like groups, rings etc. there exist maps (homomorphisms) which preserve this structure. In our specific case, a homomorphism :math:`\phi: V_{n,r} \to V_{n,s}` is a function satisfying

.. math:: w \alpha_i \phi = w \phi \alpha_i \qquad \text{and} \qquad
          w_1 \dots w_n \lambda \phi = w_1 \phi \dots w_n \phi \lambda.

As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. 

This module works with automorphisms of :math:`V_{n,r}(\boldsymbol{x})`. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.
 
.. testsetup::
	
	from thompson.automorphism import *
	from thompson.generators import Generators
	from thompson.word import Word
"""

__all__ = ["Automorphism", "Orbit"]

from collections import deque
from enum import Enum
from itertools import chain
from io import StringIO

from .word import *
from .generators import Generators
from .full_tree import FullTree

class Orbit(Enum):
	r"""Let :math:`y` be a word, and let :math:`X` be an expansion of the standard basis :math:`\boldsymbol{x}=\{x_1, \dotsc, x_n\}`; finally let :math:`\phi` be an  :class:`Automorphism`.
	
	We call the set :math:`\{ y \phi^i\}_{i\in \mathbb Z}` the :math:`\phi`-orbit of :math:`y`. In [Hig]_ (section 9), Higman showed that these :math:`\phi`-orbits could only intersect with :math:`X\langle A \rangle` in certain ways. Here, :math:`X\langle A \rangle` is the set of :meth:`simple <thompson.word.Word.is_simple>` words starting with an :math:`x \in X`.
	
	Let us refer to the intersection of the :math:`\phi`-orbit of :math:`y` as just the *orbit* of :math:`y`. This class is an :mod:`enumeration <py3:enum>` type which labels Higman's five orbit types; see section 4.1 of the paper.
	
	 1. *Complete infinite.* The orbit exists for all :math:`i \in \mathbb Z` and each element of the orbit is different.
	 2. *Complete finite.* The orbit exists for all :math:`i \in \mathbb Z`, but eventually it repeats itself.
	 3. *Right semi-infinite.* The forward orbit :math:`\{ y \phi^n\}_{n\in \mathbb N}` exists and does not repeat itself; however no part of the backward orbit :math:`\{ y \phi^{-n}\}_{n\in \mathbb N}` exists.
	 4. *Left semi-infinite.* The backward orbit exists and does not repeat itself; but no part of the forward orbit exists.
	 5. *Incomplete*. Only a finite number of images :math:`y\phi^{-n}, \dotsc, y\phi^{-1}, y, y\phi, \dotsc, y\phi^m` exist and all others do not.
	"""
	complete_infinite = 1
	complete_finite = 2
	right_semi_infinite = 3
	left_semi_infinite = 4
	incomplete = 5

def _concat(words):
	"""Takes an iterable *words* which yields lists of integers representing words. Returns a tuple containing all the *words* concatenated together, with a zero (lambda) added on the end."""
	return tuple(chain.from_iterable(words)) + (0,)

#TODO. Check the assumption that the bases consist of simple words only (no lambdas)
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
		#TODO: is this really 100% necessary? Seem to remember the big lemma/algorithm for QNF needed it.
		#How do you know that you've found the smallest possible nice basis if you haven't kept everything as small as possible throughout?
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
	
	@staticmethod
	def _reduce(domain, range):
		"""Contracts the domain generators whenever the corresponding contraction in range is possible. (This corresponds to reducing a tree pair diagram.)
			
			>>> from thompson.examples import cyclic_order_six as cyclic
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
		"""
		#TODO finish me and test me to high heaven.
		basis = self._minimal_expansion()
		print(basis)
		#expand basis until each no element's orbit has finite intersection with X<A>
		i = 0
		while i < len(basis):
			type = self._orbit_type(basis[i], basis)
			print("Orbit type:", type, '\n')
			if type is Orbit.incomplete:
				print('expand basis')
				basis.expand(i)
			else:
				i += 1
		return basis
	
	def _minimal_expansion(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \cup W`. See example 4.25.
		
		>>> from thompson.examples import *
		>>> cyclic_order_six._minimal_expansion()
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2 a1', 'x1 a1 a2 a2', 'x1 a2'])
		>>> example_4_11._minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2'])
		>>> example_4_12._minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2'])
		>>> example_4_25._minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
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
		
		>>> from thompson.examples import example_4_25 as ex
		>>> basis = ex._minimal_expansion()
		>>> ex._orbit_type(Word("x a1", 2, 1), basis)
		<Orbit.right_semi_infinite: 3>
		"""
		#TODO. Tests from the examples
		print('Forward Orbit for', y)
		right_infinite = self._test_semi_infinite(y, basis, forward=True)
		if isinstance(right_infinite, Orbit):
			return right_infinite #periodic
		print('right_infinite:', right_infinite)
		
		print('Backward orbit for', y)
		left_infinite = self._test_semi_infinite(y, basis, forward=False)
		assert not isinstance(left_infinite, Orbit), "Orbit is not periodic going forward but is going backward."
		
		if right_infinite and left_infinite:
			return Orbit.complete_infinite
		elif right_infinite and not left_infinite:
			return Orbit.right_semi_infinite
		elif not right_infinite and left_infinite:
			return Orbit.left_semi_infinite
		else:
			return Orbit.incomplete
	
	def _test_semi_infinite(self, y, basis, forward=True):
		"""Computes the orbit type of *y* with respect to *basis* in the forward direction. (Use ``forward=False`` to go backwards."""
		i = 0
		image = y
		images = [y]
		while True:
			#Compute the image y\phi^i as y\phi^{i-1} \phi
			image = self.image(image, inverse=not forward)
			print('power #', len(images), 'is', image)
			#1. Is this image in X<A>?
			if not basis.is_above(image): #not in X<A>
				return False #NOT semi_infinite in the given direction
			
			#2. Look for basis elements which are prefixes of the new image
			prefixes = [gen for gen in basis if gen.is_above(image)]
			
			#3. For each previous image:
			for previous in images:
				#Is this the same as the word we've just computed?
				if previous == image:
					return Orbit.complete_finite #Perodic and infinite in both directions.
				#Otherwise, is there a generator which is an initial segment of both the previous and current images? 
				for generator in prefixes:
					tail = generator.test_above(image)
					if tail is not None:
						#We've found a match: both *image* and *previous* start with *generator*
						return True #IS semi_infinite in the given direction
			images.append(y)


#TODO. Compose and invert automorphisms. Basically move all the functionality from TreePair over to this class and ditch trree pair.
#TODO method to decide if the automorphism is in (the equivalent of) F, T, or V.
#TODO the named elements A, B, C, X_n of Thompson's V.

