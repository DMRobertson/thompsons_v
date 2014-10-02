"""
.. testsetup::
	
	from thompson.automorphism import *
	from thompson.generators import Generators
	from thompson.word import Word
"""

__all__ = ["Automorphism"]

from collections import deque
from itertools import chain
from io import StringIO

from .word import *
from .generators import Generators
from .full_tree import FullTree

def _concat(words):
	"""Takes an iterable *words* which yields lists of integers representing words. Returns a tuple containing all the *words* concatenated together, with a zero (lambda) added on the end."""
	return tuple(chain.from_iterable(words)) + (0,)

#TODO. Check the assumption that the bases consist of simple words only (no lambdas)
class Automorphism:
	r"""Represents an automorphism of :math:`V_{n,r}` by specifying two bases. This class keeps track of the mapping between bases.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	:ivar domain: a :class:`generating set <Generators>` of preimages 
	:ivar range: a :class:`generating set <Generators>` of images.
	"""
	
	def __init__(self, arity, alphabet_size, domain, range):
		r"""Creates an automorphism, given the *arity* :math:`n` and *alphabet_size* :math:`r`. Two bases *domain* and *range* are given. The automorphism maps elements so that order is preserved:
		
			.. math:: \text{domain}_i \mapsto \text{range}_i
		
		After creation, the automorphism is reduced. This means that collections of rules like
		
			.. math:: u\alpha_1 \mapsto v\alpha_1, \dotsc, u\alpha_n \mapsto v\alpha_n
		
		are replaced by simpler rules :math:`u \mapsto v` whenever possible.
		
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
		#TODO: is this really 100% neccesary?
		Automorphism._reduce(domain, range)
		
		self.arity = arity
		self.alphabet_size = alphabet_size
		self.domain = domain
		self.range = range
		self._dict = {}
		for d, r in zip(self.domain, self.range):
			self._dict[d] = r
		#Compute and cache the images of any element in X<A> above self.domain
		for root in Generators.standard_basis(self.arity, self.alphabet_size):
			self._image_simple_above_domain(root)
	
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
	
	#Finding images of words
	def __getitem__(self, key):
		"""Computes the image of a *key* under the given automorphism. The result is cached for further usage later.
		
		The word must be given as one of:
		
		- a list of integers (see the :mod:`~thompson.word` module), or
		- a :class:`~thompson.word.Word` instance.
		
		The input need not be in standard form. This method
		
		1. Checks if the image of *key* has been cached, and returns the image if so.
		2. Checks if the image of the standard form of *key* has been cached, and returns the image if so.
		3. Computes the image of *key* under the automorphism, then caches and returns the result.
		
			>>> from thompson.examples import example_4_25
			>>> #An element of the domain---just a lookup
			>>> print(example_4_25['x1 a1'])
			x1 a1 a1 a1
			>>> #A word below a the domain words
			>>> print(example_4_25['x1 a1 a2 a2'])
			x1 a1 a1 a1 a2 a2
			>>> #Above domain words---have to expand.
			>>> print(example_4_25['x1'])
			x1 a1 a1 a1 x1 a1 a1 a2 x1 a2 a2 x1 a1 a2 L x1 a2 a1 L L L
			>>> #Let's try some words not in standard form
			>>> print(example_4_25['x1 a1 a1 x1 a1 a2 L'])
			x1 a1 a1 a1
			>>> print(example_4_25['x1 a1 a1 x1 a1 a2 L a2 a1'])
			x1 a1 a1 a1 a2 a1
			>>> print(example_4_25['x1 a1 a1 x1 a2 a2 L'])
			x1 a1 a1 a1 a1 x1 a2 a2 x1 a1 a2 L x1 a2 a1 L L
		
		:rtype: a Word instance (which are always in standard form).
		"""
		try:
			return self._dict[key]
		except KeyError:
			if isinstance(key, Word):
				word = key
			elif isinstance(key, (str, tuple)):
				word = Word(key, self.arity, self.alphabet_size)
			image = self._get_image(word)
			#TODO. Is this the right decision?
			# self._dict[key] = image
			return image
	
	def _get_image(self, word):
		"""This is the meat and bones which computes or looks up the images of a *word* in standard form. This method sets ``self._dict[word] = image``, where *image* is the image of *word* and returns *image*. 
		
		If *word* is came from tuple of *letters* not in standard form, the caller of this method may then set ``self._dict[letters] = image`` if it wishes."""
		try:
			return self._dict[word]
		except KeyError:
			pass
		
		#1. First deal with the easy words (no lambdas).
		if word.is_simple():
			#During initialisation we compute the image of everything above the domain.
			#Thus *word* must be below the domain, so we remove alphas until we find something we know the image of.
			return self._image_simple_below_domain(word)
		#2. Words with lambdas in them are going to be a bit trickier.
		return self._image_of_lambda(word)
	
	def _image_simple_above_domain(self, word):
		r"""Finds the image of a *word* in standard form above *self.domain* by expansion.
		
		.. math:: w\phi = w\alpha\lambda\phi &= w\alpha_1     \dots w\alpha_n \lambda \phi \\
		                                     &= w\alpha_1\phi \dots w\alpha_n \phi \lambda
		
		Images are cached once computed.
		"""
		try: 
			return self._dict[word]
		except KeyError:
			img_letters = _concat(self._image_simple_above_domain(child) for child in word.expand())
			#NOT in standard form.
			img_letters = standardise(img_letters, self.arity)
			image = Word(img_letters, self.arity, self.alphabet_size, preprocess=False)
			self._dict[word] = image
			return image
	
	def _image_simple_below_domain(self, word):
		r"""This method takes a :meth:`simple <Word.is_simple>` *word* of the form :math:`d \alpha_i_1 \dots \alpha_i_m` where :math:`d` is the largest such word whose image has already been computed. The images of the words
		
		.. math:: d \alpha_i_1, d \alpha_i_1 \alpha_i_2, \dotsc, d\alpha_i_1 \alpha i_2 \dots \alpha i_m
		
		are then computed and cached. The final image in this list (i.e. that of the original *word* is returned).
		"""
		i = 1
		while True:
			head, tail = word.split(i)
			if head in self._dict:
				break
			i += 1
		head = Word(head, self.arity, self.alphabet_size, preprocess=False)
		image = self._dict[head] #is a word
		for _ in range(i):
			alpha, tail = tail[0], tail[1:]
			head = head.alpha(-alpha) #.alpha() returns a word
			image = image.alpha(-alpha) #.alpha() returns a word
			self._dict[head] = image
		assert len(tail) == 0
		return image
		
	def _image_of_lambda(self, letters, in_standard_form=True):
		r"""Let *letters* be a list of integers describing a word which ends in a lambda. If *in_standard_form* is True, then the *letters* should describe a word in standard form. Otherwise, with the exception that words of the form
		
		.. math:: w\alpha\lambda = w \alpha_1 w \alpha_2 \dots w\alpha_n \lambda
		
		are permitted.
		
		This method expands the word into the arguments of the lambda, then calls :meth:`_get_image` on each subword; the process continues recursively. Once we have computed/retreived all of the images we need, the images are concatenated and standardised, forming the image of *letters*.
		
		If *letters* is *in_standard_form*, then the image is cached before returning. If not, it is the responsibility of the caller to set ``self._dict[word] = image``, where *word* is the standardised form of *letters*.
		
		:raises ValueError: if the last letter in *word* is not a lambda.
		"""
		subwords = lambda_arguments(letters)
		letters = _concat(self._get_image(word) for word in subwords)
		letters = Word(letters, self.arity, self.alphabet_size)
		image = Word(letters, self.arity, self.alphabet_size, preprocess=False)
		if in_standard_form:
			self._dict[letters] = image
		return image
	
	def minimal_expansion(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \cup W`. See example 4.25.
		
		>>> from thompson.examples import cyclic_order_six, example_4_25
		>>> cyclic_order_six.minimal_expansion()
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2 a1', 'x1 a1 a2 a2', 'x1 a2'])
		>>> example_4_25.minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
		"""
		basis = Generators.standard_basis(self.arity, self.alphabet_size)
		i = 0
		while i < len(basis):
			b = basis[i]
			if b in self.domain or b in self.range:
				i += 1
			else:
				basis.expand(i)
		return basis
	
	def to_quasinormal_form(self):
		"""The plan.
		
		1. Reduce the automorphism (eliminate carets) - DONE 
		2. Find elements above (Y union W) - DONE indirectly, see 3.
		3. Expand std basis until it contains none of the elements from step 2. DONE
		4. Test each element of this expanded basis to see if the have nice orbits.
			a. If an element doesn't, expand it and try again with its children.
		5. When everything in the basis has nice orbits we are done.
		"""
		pass
	
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
			value = self[key]
			output.write(fmt.format(key, value))
		return output.getvalue()
	
	def dump_mapping(self, **kwargs):
		print('Automorphism of V_{},{} specified by {} generators (after reduction).'.format(
		  self.arity, self.alphabet_size, len(self.domain)), **kwargs)
		print('The generators are marked by an asterisk.', **kwargs)
		print('Following from the above {} rules, the following computations have been cached.'.format(
		  len(self.domain)), **kwargs)
		max_len = 0
		for key in self._dict:
			max_len = max(max_len, len(str(key)))
		fmt = "{!s: <" + str(max_len) + "} -> {!s}"
		for key in sorted(self._dict):
			prefix = " * " if key in self.domain else "   "
			value = self[key]
			print(prefix + fmt.format(key, value), **kwargs)
		

#TODO. Compose and invert automorphisms. Basically move all the functionality from TreePair over to this class and ditch trree pair.

