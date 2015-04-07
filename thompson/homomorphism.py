r"""The algebra of :mod:`words <thompson.word>` has its own structure, and just like groups, rings etc. there exist homomorphisms which preserve this structure. In our specific case, a homomorphism :math:`\phi: V_{n,r} \to V_{n,s}` is a function satisfying

.. math:: w \alpha_i \phi = w \phi \alpha_i \qquad \text{and} \qquad
          w_1 \dots w_n \lambda \phi = aw_1 \phi \dots w_n \phi \lambda.

In other words, a homomorphism is a function which 'commutes' with the algebra operations :math:`\alpha_i` and :math:`\lambda`.

.. testsetup::
    
    from thompson.generators   import Generators
    from thompson.homomorphism import *
    from thompson.automorphism import Automorphism
    from thompson.examples     import *
"""

from io import StringIO
from itertools import chain
import re
import sys

from .word import *
from .generators import Generators
from copy import copy

__all__ = ['Homomorphism']

class Homomorphism:
	r"""Let :math:`f: D \to R` be some map embedding a basis :math:`D` for :math:`V_{n,r}` into another algebra :math:`V_{n,s}` of the same :class:`~thompson.word.Signature`. This map uniquely extends to a homomorphism of algebras :math:`\psi : V_{n,r} \to V_{n,s}`.
	
	:ivar domain: a :class:`basis <thompson.generators.Generators>` of preimages 
	:ivar range:  a :class:`basis <thompson.generators.Generators>` of images.
	
	The next few attributes are used internally when constructing free factors. We need to have a way to map back to the parent automorphisms.
	
	:ivar domain_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	:ivar range_relabeller: the same.
	"""
	
	#Initialisation
	def __init__(self, domain, range, reduce=True):
		"""Creates a homomorphism with the specified *domain* and *range*. Sanity checks are performed so that the arguments do genuinely define a basis.
		
		By default, any redudancy in the mapping is reduced. For example, the rules ``x1 a1 a1 -> x1 a2 a1`` and ``x1 a1 a2 -> x1 a2 a2`` reduce to the single rule ``x1 a1 -> x1 a2``. The *reduce* argument can be used to disable this. This option is intended for internal use only. [#footnote_why_optional_reduce]_
		
		:raises TypeError: if the bases are of different sizes.
		:raises TypeError: if the algebras described by *domain* and *range* have different arities.
		:raises ValueError: if *domain* is :meth:`not a basis <thompson.generators.Generators.is_basis>`.
		"""
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
		self.domain_relabeller = None
		self.range_relabeller = None
		
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
			>>> print(load_example('cyclic_order_six'))
			PeriodicAut: V(2, 1) -> V(2, 1) specified by 5 generators (after expansion and reduction).
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
	
	#Computing Images
	def _set_image(self, d, r, sig_in, sig_out, cache):
		"""Stores the rule that phi(d) = r in the mapping dictionary and ensure that d and r are both Words."""
		assert d in sig_in, repr(d)
		assert r in sig_out, repr(r)
		cache[d] = r
	
	#Input/Output
	@classmethod
	def from_file(cls, filename):
		"""Reads in a file specifying a homomorphism and returns the homomorphism being described. Here is an example of the format::
		
			5
			(2,1)		->	(2,1)
			x1 a1 a1 a1	->	x1 a1 a1
			x1 a1 a1 a2	->	x1 a1 a2 a1
			x1 a1 a2	->	x1 a1 a2 a2
			x1 a2 a1	->	x1 a2 a2
			x1 a2 a2	->	x1 a2 a1
			This is an example intended to illustrate the format for reading and writing automorphisms to disk.
		
		There are three components
		
		- number :math:`k` of generators in domain and range
		- signatures of domain and range
		- list of :math:`k` rules domain -> range
		
		Any lines after this are ignored, and can be treated as comments. Comments are read in and added to the __doc__ attribute of the homomorphism that gets created.
		
		.. seealso: The thompson/examples directory in the source code.
		"""
		with open(filename, encoding='utf-8') as f:
			num_generators = int(f.readline())
			params = extract_signatures.match(f.readline().strip()).groups()
			d = Generators([int(params[0]), int(params[1])])
			r = Generators([int(params[2]), int(params[3])])
			for i in range(num_generators):
				d_word, r_word = (word.strip() for word in f.readline().split('->'))
				try:
					d.append(d_word)
					r.append(r_word)
				except Exception as e:
					extra = "Problem reading rule {} from line {} of {}. The original error was:\n\t".format(
					  i+1, i+3, f.name)
					raise type(e)(extra + str(e)).with_traceback(sys.exc_info()[2])
			hom = cls(d, r)
			hom.__doc__ = f.read()
			if len(hom.__doc__.strip()) == 0:
				del hom.__doc__
		return hom
	
	def save_to_file(self, filename=None):
		"""Takes a homomorphism and saves it to the file with the given *filename*. The homomorphism is stored in a format which is compatible with :meth:`from_file`.
		
		.. doctest::
			
			>>> before = random_automorphism()
			>>> before.save_to_file('test_saving_homomorphism.aut')
			>>> after = Automorphism.from_file('test_saving_homomorphism.aut')
			>>> before == after, before is after
			(True, False)
		
		.. testcleanup::
			
			import os
			os.remove('test_saving_homomorphism.aut')
			
		"""
		#todo docstring
		#todo test that this is an inverse to from_file
		if filename is None:
			clsname = type(self).__name__ 
			filename = clsname + '.' + clsname.lower()[:3]
		with open(filename, 'wt', encoding='utf-8') as f:
			print(len(self.domain), file=f)
			print("{} -> {}".format(self.domain.signature, self.range.signature), file=f)
			for d, r in zip(self.domain, self.range):
				print("{} -> {}".format(d, r), file=f)
	
	#Simple operations on homomorphisms
	def __eq__(self, other):
		"""We can test for equality of homomorphisms by using Python's equality operator ``==``.
		
			>>> phi = random_automorphism()
			>>> phi == phi
			True
			>>> phi * ~phi == Homomorphism.identity(phi.signature)
			True
		"""
		return all(self.image(w) == other.image(w) for w in chain(self.domain, other.domain))
	
	def __mul__(self, other): #self * other is used for the (backwards) composition self then other
		r"""If the current automorphism is :math:`\psi` and the *other* is :math:`\phi`, multiplication forms the composition :math:`\psi\phi`, which maps :math:`x \mapsto x\psi \mapsto (x\psi)\phi`.
		
		:raises TypeError: if the homomorphisms cannot be composed in the given order; i.e. if ``self.range.signature != other.domain.signature``.
		
		:rtype: an :class:`~thompson.mixed.MixedAut` if possible; otherwise a :class:`Homomorphism`.
		
		.. doctest::
			
			>>> phi = load_example('alphabet_size_two')
			>>> print(phi * phi)
			MixedAut: V(3, 2) -> V(3, 2) specified by 8 generators (after expansion and reduction).
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
		
		range = other.image_of_set(self.range)
		domain = copy(self.domain)
		from .automorphism import Automorphism
		for cls in (Automorphism, Homomorphism):
			try:
				return cls(domain, range)
			except Exception as e:
				print(e)
		raise e
	
	@classmethod
	def identity(cls, signature):
		"""Creates a new homo/automorphism which is the identity map on the algebra with the given *signature*.
		
			>>> print(Homomorphism.identity((3, 2)))
			PeriodicAut: V(3, 2) -> V(3, 2) specified by 2 generators (after expansion and reduction).
			x1 -> x1
			x2 -> x2
		"""
		d = Generators.standard_basis(signature)
		r = Generators.standard_basis(signature)
		from .periodic import PeriodicAut
		return PeriodicAut(d, r)
	
	def is_identity(self):
		"""Returns True if this automorphism is the identity map on the algebra with the given *signature*. Otherwise returns False.
		
			>>> aut = Homomorphism.identity(random_signature())
			>>> aut.is_identity()
			True
			>>> load_example('example_5_15').is_identity()
			False
		"""
		basis = Generators.standard_basis(self.signature)
		return self.image_of_set(basis) == basis
	
	#Finding images of words
	def image(self, key, sig_in=None, sig_out=None, cache=None):
		"""Computes the image of a *key* under the given homomorphism. The result is cached for further usage later.
		
		The *key* must be given as one of:
		
		- a string to be passed to :func:`~thompson.word.from_string`;
		- a sequence of integers (see the :mod:`~thompson.word` module); or
		- a :class:`~thompson.word.Word` instance.
		
		.. note::
			
			All other arguments are intended for internal use only. [#footnote_why_optional_image_args]_
		
		The input need not be in standard form. This method
		
		1. Converts *key* to standard form if necessary.
		2. Checks if the image of the standard form of *key* has been cached, and returns the image if so.
		3. Computes the image of *key* under the homomorphism, then caches and returns the result.
		
		.. doctest::
		
			>>> #An element of the domain---just a lookup
			>>> phi = load_example('example_5_15')
			>>> print(phi.image('x1 a1'))
			x1 a1 a1 a1
			>>> #A word below a the domain words
			>>> print(phi.image('x1 a1 a2 a2'))
			x1 a1 a1 a1 a2 a2
			>>> #Above domain words---have to expand.
			>>> print(phi.image('x1'))
			x1 a1 a1 a1 x1 a1 a1 a2 x1 a2 a2 x1 a1 a2 L x1 a2 a1 L L L
			>>> #Let's try some words not in standard form
			>>> print(phi.image('x1 a1 a1 x1 a1 a2 L'))
			x1 a1 a1 a1
			>>> print(phi.image('x1 a1 a1 x1 a1 a2 L a2 a1'))
			x1 a1 a1 a1 a2 a1
			>>> print(phi.image('x1 a1 a1 x1 a2 a2 L'))
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
		"""Computes the image of a list of :class:`~thompson.generators.Generators` under the current homomorphism. The order of words in the list is preserved.
		
		:rtype: another list of :class:`~thompson.generators.Generators`.
		
		.. doctest::
			
			>>> basis = Generators.standard_basis((2,1))
			>>> basis.expand_to_size(8); print(basis)
			[x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2 a1, x1 a1 a2 a2, x1 a2 a1 a1, x1 a2 a1 a2, x1 a2 a2 a1, x1 a2 a2 a2]
			>>> print(load_example('example_5_3').image_of_set(basis))
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
		
			>>> print(load_example('cyclic_order_six'))
			PeriodicAut: V(2, 1) -> V(2, 1) specified by 5 generators (after expansion and reduction).
			x1 a1 a1    -> x1 a1 a1      
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2         
			x1 a2 a1    -> x1 a1 a2 a2 a1
			x1 a2 a2    -> x1 a1 a2 a1   
		"""
		relabelled = self.domain_relabeller is not None
		output = StringIO()
		output.write(self._string_header())
		
		if relabelled:
			output.write("\nThis automorphism was derived from a parent automorphism.\n'x' and 'y' represent root words of the parent and current derived algebra, respectively.")
			domain_relabelled, range_relabelled = self.relabel()
			rows = format_table(
				domain_relabelled, self.domain, self.range, range_relabelled,
				sep = ['~>   ', '=>', '   ~>'], root_names = 'xyyx')
		else:
			rows = format_table(self.domain, self.range)
		
		for row in rows:
			output.write('\n')
			output.write(row)
		return output.getvalue()
	
	#Relabelling
	def add_relabellers(self, domain_relabeller, range_relabeller):
		""":raises LookupError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range*
		:raises TypeError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range* (as appropriate)
		"""
		if (domain_relabeller is None) != (range_relabeller is None):
			raise LookupError("Must either specify both relabellers or neither.", domain_relabeller, range_relabeller)
		
		if domain_relabeller is not None and domain_relabeller.domain.signature != self.domain.signature:
			raise TypeError('Domain relabeller signature {} does not match automorphism\'s domain signature {}.'.format(
			  domain_relabeller.domain.signature, domain.signature))
		
		if range_relabeller is not None and range_relabeller.domain.signature != self.range.signature:
			raise TypeError('Range relabeller signature {} does not match automorphism\'s range signature {}.'.format(
			  range_relabeller.domain.signature, range.signature))
		
		self.domain_relabeller = domain_relabeller
		self.range_relabeller = range_relabeller
	
	def relabel(self):
		r"""If this automorphism was derived from a parent automorphism, this converts back to the parent algebra after doing computations in the derived algebra.
		
		In the following example :meth:`~thompson.mixed.MixedAut.test_conjugate_to` takes a pure periodic automorphism and extracts factors. A conjugator :math:`\rho` is produced by :meth:`the overridden version of this method <thompson.periodic.PeriodicAut.test_conjugate_to>`. Finally :math:`\rho` is relabelled back to the parent algebra.
		
		:raises AttributeError: if the factor has not been assigned any relabellers.
		
		.. doctest::
			
			>>> psi, phi = load_example_pair('example_5_12')
			>>> rho = psi.test_conjugate_to(phi)
			>>> print(rho)
			PeriodicAut: V(2, 1) -> V(2, 1) specified by 6 generators (after expansion and reduction).
			x1 a1 a1 a1 a1 -> x1 a1 a2   
			x1 a1 a1 a1 a2 -> x1 a2 a2   
			x1 a1 a1 a2    -> x1 a1 a1 a1
			x1 a1 a2       -> x1 a2 a1 a1
			x1 a2 a1       -> x1 a1 a1 a2
			x1 a2 a2       -> x1 a2 a1 a2
			>>> rho * phi == psi * rho
			True
		"""
		if self.domain_relabeller is None or self.range_relabeller is None:
			raise AttributeError("This factor has not been assigned relabellers.")
		return self.domain_relabeller.image_of_set(self.domain), self.range_relabeller.image_of_set(self.range)

def format_table(*columns, sep=None, root_names=None):
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