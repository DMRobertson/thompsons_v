r"""The algebra of :mod:`words <thompson.word>` has its own structure, and just like groups, rings etc. there exist homomorphisms which preserve this structure. In our specific case, a homomorphism :math:`\phi: V_{n,r} \to V_{n,s}` is a function satisfying

.. math:: w \alpha_i \phi = w \phi \alpha_i \qquad \text{and} \qquad
          w_1 \dots w_n \lambda \phi = w_1 \phi \dots w_n \phi \lambda.

In other words, a homomorphism is a function which 'commutes' with the algebra operations :math:`\alpha_i` and :math:`\lambda`.

.. testsetup::
    
    from thompson.homomorphism import *
    from thompson.examples import *
"""

from io import StringIO
from itertools import chain
import re

from .word import *
from .generators import Generators
from copy import copy

#Extracted the bits responsible for defining a homomorphism from the automorphism class.
class Homomorphism:
	r"""Let :math:`f: D \to R` be some map embedding a basis :math:`D` for :math:`V_{n,r}` into another algebra :math:`V_{n,s}` of the same :class:`~thompson.word.Signature`. This map uniquely extends to a homomorphism of algebras :math:`\psi : V_{n,r} \to V_{n,s}`.
	
	:ivar domain: a :class:`basis <thompson.generators.Generators>` of preimages 
	:ivar range:  a :class:`set <thompson.generators.Generators>` of images.
	
	:raises TypeError: if the bases are of different sizes.
	:raises TypeError: if the algebras described by *domain* and *range* have different arities.
	:raises ValueError: if *domain* is :meth:`not a basis <thompson.generators.Generators.is_basis>`.
	"""
	#TODO docstring for reduce
	def __init__(self, domain, range, reduce=True):
		#The boring checks
		if len(domain) != len(range):
			raise TypeError("Domain basis has {} elements, but range basis has {} elements.".format(
			  len(domain), len(range)))
		
		if domain.signature.arity != range.signature.arity:
			raise TypeError("Domain arity {} does not equal the range arity {}.".format(
			  domain.signature.arity, range.signature.arity))
		
		#Expand any non-simple words
		Homomorphism._expand(domain, range)
		domain, range = Generators.sort_mapping_pair(domain, range)
		
		#Remove any redundancy---like reducing tree pairs.
		#How do you know that you've found the smallest possible nice basis if you haven't kept everything as small as possible throughout?
		if reduce:
			Homomorphism._reduce(domain, range)
		
		#Check that domain is a free generating set
		i, j = domain.test_free()
		if not(i == j == -1):
			raise ValueError("Domain is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		#Check to see that the domain generates all of V_{n,r}
		missing = domain.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Domain {} does not generate V_{}. Missing elements are {}.".format(
			  domain, domain.signature, [format(x) for x in missing]))
		
		self.domain = domain
		self.range = range
		
		#Setup the mapping cache
		self._map = {}
		for d, r in zip(self.domain, self.range):
			self._set_image(d, r, self.domain.signature, self.range.signature, self._map)
			
		#Compute and cache the images of any simple word above self.domain.
		for root in Generators.standard_basis(domain.signature):
			self._image_simple_above_domain(root, self.domain.signature, self.range.signature, self._map)
	
	@staticmethod
	def _expand(domain, range):
		r"""Expands the pair of generating sets where necessary until all words in both sets are simple.
		
			>>> g = Generators((2, 2), ['x1 a1 a1', 'x1 a2 x1 a1 L', 'x2 x1 L', 'x2 a2']);
			>>> h = Generators((2, 2), ['x1 a1 x2 L', 'x2 a1 a1 a2', 'x2 x1 L', 'x2 a2']);
			>>> Homomorphism._expand(g, h)
			>>> print(g, h, sep='\n')
			[x1 a1 a1 a1, x1 a1 a1 a2, x1 a2, x1 a1, x2, x1, x2 a2]
			[x1 a1, x2, x2 a1 a1 a2 a1, x2 a1 a1 a2 a2, x2, x1, x2 a2]
		"""
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
			
			>>> #This is given by 6 generators, but after reduction consists of 5:
			>>> print(cyclic_order_six)
			Automorphism: V(2, 1) -> V(2, 1) specified by 5 generators (after expansion and reduction).
			x1 a1 a1    -> x1 a1 a1      
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2         
			x1 a2 a1    -> x1 a1 a2 a2 a1
			x1 a2 a2    -> x1 a1 a2 a1   
			>>> #Swaps x1 and x2. Should reduce to two generators instead of 3
			>>> domain = Generators((2, 2), ["x1", "x2 a1", "x2 a2"])
			>>> range  = Generators((2, 2), ["x2", "x1 a1", "x1 a2"])
			>>> print(Homomorphism(domain, range))
			Homomorphism: V(2, 2) -> V(2, 2) specified by 2 generators (after expansion and reduction).
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
	
	def _set_image(self, d, r, sig_in, sig_out, cache):
		"""Stores the rule that phi(d) = r in the mapping dictionary and ensure that d and r are both Words."""
		assert d in sig_in, repr(d)
		assert r in sig_out, repr(r)
		cache[d] = r
	
	#Alternative constructors
	@classmethod
	def from_file(cls, filename):
		"""Reads in a file specifying an automorphism and returns a pair *(aut, name)*. Here is an example of the format::
		5
		(2,1)		->	(2,1)
		x1 a1 a1 a1	->	x1 a1 a1
		x1 a1 a1 a2	->	x1 a1 a2 a1
		x1 a1 a2	->	x1 a1 a2 a2
		x1 a2 a1	->	x1 a2 a2
		x1 a2 a2	->	x1 a2 a1
		
		- number of generators
		- signatures of domain and range
		- list of rules domain -> range
		"""
		with open(filename, encoding='utf-8') as f:
			num_generators = int(f.readline())
			params = extract_signatures.match(f.readline().strip()).groups()
			d = Generators([int(params[0]), int(params[1])])
			r = Generators([int(params[2]), int(params[3])])
			for _ in range(num_generators):
				d_word, r_word = (word.strip() for word in f.readline().split('->'))
				d.append(d_word)
				r.append(r_word)
		return cls(d, r)
	
	#Simple operations on homomorphisms
	def __eq__(self, other):
		return all(self.image(w) == other.image(w) for w in chain(self.domain, other.domain))
	
	def __mul__(self, other): #self * other is used for the (backwards) composition self then other
		r"""If the current automorphism is :math:`\psi` and the *other* is :math:`\phi`, multiplication forms the composition :math:`\psi\phi`, which maps :math:`x \mapsto x\psi \mapsto (x\psi)\phi`.
		
		:raises TypeError: if the homomorphisms cannot be composed in the given order; i.e. if ``self.range.signature != other.domain.signature``.
		
		:rtype: an :class:`~thompson.automorphism.Automorphism` if possible; otherwise a :class:`Homomorphism`.
		
		.. doctest::
			
			>>> print(alphabet_size_two * alphabet_size_two)
			Automorphism: V(3, 2) -> V(3, 2) specified by 8 generators (after expansion and reduction).
			x1 a1    -> x1 a1      
			x1 a2    -> x1 a2 a3 a3
			x1 a3 a1 -> x1 a3      
			x1 a3 a2 -> x1 a2 a2   
			x1 a3 a3 -> x1 a2 a1   
			x2 a1    -> x2         
			x2 a2    -> x1 a2 a3 a2
			x2 a3    -> x1 a2 a3 a1
		"""
		if not isinstance(other, Homomorphism):
			return NotImplemented
		if self.range.signature != other.domain.signature:
			raise TypeError("Signatures {}->{} and {}->{} do not match.".format(
			  self.domain.signature, self.range.signature, other.domain.signature, other.range.signature))
		
		range = Generators(other.range.signature)
		range = other.image_of_set(self.range)
		domain = copy(self.domain)
		
		from thompson.automorphism import Automorphism
		type = Automorphism if self.domain.signature == other.range.signature else Homomorphism
		return type(domain, range)
	
	#Finding images of words
	def image(self, key, sig_in=None, sig_out=None, cache=None):
		"""Computes the image of a *key* under the given homomorphism. The result is cached for further usage later.
		
		The *key* must be given as one of:
		
		- a sequence of integers (see the :mod:`~thompson.word` module), or
		- a :class:`~thompson.word.Word` instance.
		
		The input need not be in standard form. This method
		
		1. Converts *key* to standard form if necessary.
		2. Checks if the image of the standard form of *key* has been cached, and returns the image if so.
		3. Computes the image of *key* under the homomorphism, then caches and returns the result.
		
		.. doctest::
		
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
		
		:rtype: a :class:`~thompson.word.Word` instance (which are always in standard form).
		"""
		sig_in = sig_in or self.domain.signature
		sig_out = sig_out or self.range.signature
		cache = cache or self._map
		try:
			return cache[key]
		except KeyError:
			return self._compute_image(key, sig_in, sig_out, cache)
	
	def _compute_image(self, key, sig_in, sig_out, cache):
		if isinstance(key, Word):
			if key.signature != sig_in:
				raise TypeError("Signature {} of input {} does not match {}".format(
				  key.signature, key, sig_in))
			word = key
		elif isinstance(key, (str, tuple)):
			word = Word(key, sig_in)
		
		try:
			return cache[word]
		except KeyError:
			pass
		
		#1. First deal with the easy words (no lambdas).
		if word.is_simple():
			#Cache should already have dealt with anything above the domain
			return self._image_simple_below_domain(word, sig_in, sig_out, cache)
		#2. Words with lambdas in them are going to be a bit trickier.
		return self._image_of_lambda(word, sig_in, sig_out, cache)
	
	def _image_simple_above_domain(self, word, sig_in, sig_out, cache):
		r"""Finds the image of a *word* in standard form above *self.domain* by expansion.
		
		.. math:: w\phi = w\alpha\lambda\phi &= w\alpha_1     \dots w\alpha_n \lambda \phi \\
		                                     &= w\alpha_1\phi \dots w\alpha_n \phi \lambda
		
		Images are cached once computed.
		"""
		assert word.signature == sig_in
		try: 
			return cache[word]
		except KeyError:
			img_letters = _concat(
			  self._image_simple_above_domain(child, sig_in, sig_out, cache)
			    for child in word.expand())
			#NOT in standard form.
			img_letters = standardise(img_letters, sig_in)
			image = Word(img_letters, sig_out, preprocess=False)
			self._set_image(word, image, sig_in, sig_out, cache)
			return image
	
	def _image_simple_below_domain(self, word, sig_in, sig_out, cache):
		r"""This method takes a :meth:`simple <Word.is_simple>` *word* of the form :math:`d \alpha_i_1 \dots \alpha_i_m` where :math:`d` is the largest such word whose image has already been computed. The images of the words
		
		.. math:: d \alpha_i_1, d \alpha_i_1 \alpha_i_2, \dotsc, d\alpha_i_1 \alpha i_2 \dots \alpha i_m
		
		are then computed and cached. The final image in this list (i.e. that of the original *word* is returned).
		"""
		i = 1
		while True:
			head, tail = word.rsplit(i)
			if head in cache:
				break
			i += 1
		head = Word(head, sig_in, preprocess=False)
		image = cache[head] #is a word
		for _ in range(i):
			alpha, tail = tail[0], tail[1:]
			head = head.alpha(-alpha) #.alpha() returns a word
			image = image.alpha(-alpha) #.alpha() returns a word
			self._set_image(head, image, sig_in, sig_out, cache)
		assert len(tail) == 0
		return image
		
	def _image_of_lambda(self, word, sig_in, sig_out, cache):
		r"""This method accepts a *word* which ends in a lambda and extracts the arguments of the lambda. Next, the method calls :meth:`_get_image` on each argument, continuing the process recursively. Once we have computed/retreived all of the images we need, the images are concatenated and standardised, forming the image of *word*. This image is cached and returned.
		
		:raises ValueError: if the last letter in *word* is not a lambda.
		"""
		subwords = lambda_arguments(word)
		letters = _concat(self._compute_image(word, sig_in, sig_out, cache) for word in subwords)
		letters = standardise(letters, sig_out)
		image = Word(letters, sig_out, preprocess=False) #can skip validation
		self._set_image(word, image, sig_in, sig_out, cache)
		return image
	
	def image_of_set(self, set, sig_in=None, sig_out=None, cache=None):
		"""Computes the image of a set of :class:`~thompson.generators.Generators` under the current homomorphism.
		
		:rtype: another set of :class:`~thompson.generators.Generators`.
		
		.. doctest::
			
			>>> basis = Generators.standard_basis((2,1))
			>>> basis.expand_to_size(8); print(basis)
			[x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2 a1, x1 a1 a2 a2, x1 a2 a1 a1, x1 a2 a1 a2, x1 a2 a2 a1, x1 a2 a2 a2]
			>>> print(example_5_3.image_of_set(basis))
			[x1 a1 a1 a1 x1 a1 a1 a2 a1 L, x1 a1 a1 a2 a2, x1 a1 a2 a2, x1 a1 a2 a1, x1 a2 a1 a1 a1, x1 a2 a1 a1 a2, x1 a2 a1 a2, x1 a2 a2]
		"""
		sig_in = sig_in or self.domain.signature
		sig_out = sig_out or self.range.signature
		cache = cache or self._map
		
		if set.signature != sig_in:
			raise ValueError("Set signature {} does not match the input signature {}.".format(
			  set.signature, sig_in))
		images = Generators(sig_out)
		for preimage in set:
			images.append(self._compute_image(preimage, sig_in, sig_out, cache))
		
		return images
	
	#Printing
	def _string_header(self):
		return "{}: V{} -> V{} specified by {} generators (after expansion and reduction).".format(
		  type(self).__name__, self.domain.signature, self.range.signature, len(self.domain))
	
	def __str__(self):
		"""Printing an automorphism gives its arity, alphabet_size, and lists the images of its domain elements.
		
			>>> print(cyclic_order_six)
			Automorphism: V(2, 1) -> V(2, 1) specified by 5 generators (after expansion and reduction).
			x1 a1 a1    -> x1 a1 a1      
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2         
			x1 a2 a1    -> x1 a1 a2 a2 a1
			x1 a2 a2    -> x1 a1 a2 a1   
		"""
		output = StringIO()
		output.write(self._string_header())
		rows = self._format_table(self.domain, self.range)
		for row in rows:
			output.write('\n')
			output.write(row)
		return output.getvalue()
	
	@staticmethod
	def _format_table(*columns, sep=None, root_names=None):
		for row in zip(*columns):
			break
		num_columns = len(row)
		
		if sep is None:
			sep = [' -> '] * (num_columns - 1)
		else:
			assert len(sep) == (num_columns - 1)
			sep = [' {} '.format(s) for s in sep]
		
		if root_names is None:
			root_names = 'x' * num_columns
		else:
			assert len(root_names) == num_columns
		
		max_width = [0] * num_columns
		for row in zip(*columns):
			for i, entry in enumerate(row):
				max_width[i] = max(max_width[i], len(str(entry)))
		
		column = "{{!s: <{}}}"
		fmt = ""
		for i, width in enumerate(max_width):
			if i > 0:
				fmt += sep[i - 1]
			fmt += column.format(max_width[i])
		
		for row in zip(*columns):
			row = [str(entry).replace('x', root_names[i]) for i, entry in enumerate(row)]
			yield fmt.format(*row)

#Used in from_file()
extract_signatures = re.compile(r"""
	\( \s* (\d+) [\s,]+ (\d+) \s* \)
	[\s\->]+
	\( \s* (\d+) [\s,]+ (\d+) \s* \)""", re.VERBOSE)