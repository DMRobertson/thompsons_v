"""
.. testsetup:: 
	
	from thompson.word import Word, format
	from thompson.generators import *
	from thompson.examples import *
"""
from copy import copy
from collections import defaultdict

from . import word
from .word import Word, Signature

__all__ = ["Generators", "rewrite_set", "rewrite_word"]

class Generators(list):
	r"""An ordered subset of :math:`V_{n, r}`, together with methods which treat such sets as generating sets and bases. Internally, this is a subclass of :class:`list <py3:list>`, so it is possible to reorder and otherwise modify collections of Generators.
	
	:ivar signature: The :class:`~thompson.word.Signature` of the algebra this set belongs to.
	"""
	
	def __init__(self, signature, generators=None):
		"""When creating a generating set, you must specify the algebra :math:`V_{n, r}` it is a subset of. A list of generators can optionally be provided. Each generator is passed to :meth:`append`.
		"""
		if not isinstance(signature, Signature):
			signature = Signature(*signature)
		self.signature = signature
		
		self.cache = None
		
		if generators is not None:
			for g in generators:
				self.append(g)
	
	def append(self, w):
		"""Adds *w* to this generating set. If *w* is a string, a :class:`~thompson.word.Word` object is created and assigned the same *arity* and *alphabet_size* as the generating set.
		
		:raises TypeError: if *w* is neither a string nor a Word.
		:raises IndexError: if *w* has a different arity to the generating set.
		:raises IndexError: if *w* has a larger alphabet_size the generating set.
		:raises ValueError: if *w* is already contained in this generating set.
		"""
		if isinstance(w, str):
			w = Word(w, self.signature)
		elif not isinstance(w, Word):
			raise TypeError("{:r} is neither a string nor a Word.".format(w))
		
		if self.signature != w.signature:
			raise IndexError("Can't add {} with arity {} to generating set with arity {}."
			  .format(w, w.signature, self.signature))
		
		if w not in self:
			super().append(w)
		else:
			raise ValueError("The word {} is already in this set.".format(w))
	
	def __str__(self):
		return "[" + ", ".join(format(w) for w in self) + "]"
	
	def __repr__(self):
		return "Generators({}, [{}])".format(
		  self.signature, ", ".join(repr(format(w)) for w in self))
	
	def __eq__(self, other):
		if not isinstance(other, Generators):
			return NotImplemented
		return self.signature == other.signature and super().__eq__(other)
	
	#Tests on generating sets
	def test_free(self):
		r"""Tests to see if the current generating set :math:`X` is a free generating set. If the test fails, returns the first pair of indices :math:`(i, j)` found for which :math:`X_i` :meth:`is above <thompson.word.Word.is_above>` :math:`X_j`. If the test passes, returns :math:`(-1, -1)`.
		
			>>> g = Generators((2, 3), ["x1 a1", "x1 a2 a1", "x1 a2 a1 a1", "x1 a2 a2"])
			>>> g.test_free()
			(1, 2)
			>>> print(g[1], g[2], sep='\n')
			x1 a2 a1
			x1 a2 a1 a1
		
		.. seealso::  Lemma 3.16.3 of the paper.
		"""
		for i in range(len(self)):
			for j in range(i+1, len(self)):
				if self[i].is_above(self[j]):
					return i, j
				elif self[j].is_above(self[i]):
					return j, i
		return -1, -1
		
	def is_free(self):
		"""Returns True if this is a free generating set; otherwise False.
		
			>>> g = Generators((2, 3), ['x1', 'x2']);
			>>> g.is_free()
			True
			>>> g.append('x3'); g.is_free()
			True
			>>> g.append('x2 a1'); g.is_free()
			False
		"""
		return self.test_free() == (-1, -1)
	
	def simple_words_above(self):
		r"""An iterator that yields the simple words which can be obtained by contracting this basis :math:`n \geq 0` times. (So the 'above' condition is not strict.)
		
			>>> g = Generators((3, 1), ["x a1 a1 a3", "x a1 a1 a1", "x a1 a2", "x a1 a1 a2", "x a1 a3"])
			>>> for word in g.simple_words_above():
			... 	print(format(word))
			x1 a1 a1 a1
			x1 a1 a1 a2
			x1 a1 a1 a3
			x1 a1 a2
			x1 a1 a3
			x1 a1 a1
			x1 a1
		"""
		words = sorted(self) #Creates a shallow copy and then sorts it.
		yield from words
		i = 0
		arity = self.signature.arity
		while i <= len(words) - arity:
			prefix = word.are_contractible(words[i : i + arity])
			if prefix:
				yield prefix
				words[i : i + arity] = [prefix] if prefix not in words else []
				i -= (arity - 1) 
				i = max(i, 0)
			else:
				i += 1
	
	def test_generates_algebra(self):
		r"""Tests to see if this set generates all of :math:`V_{n,r}`. The generating set is sorted, and then contracted as much as possible. Then we check to see which elements of the :meth:`standard basis <standard_basis>` :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` are not included in the contraction.
		
		The test fails if there is at least one such element; in this case, the list of all such elements is returned. Otherwise the test passes, and an empty list is returned.
		
			>>> #A basis for V_2,1
			>>> Y = Generators((2, 1), ["x a1", "x a2 a1", "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
			>>> Y.test_generates_algebra()
			[]
			>>> #The same words do not generate V_2,3
			>>> Y = Generators((2, 3), ["x a1", "x a2 a1", "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
			>>> missing = Y.test_generates_algebra(); missing
			[Word('x2', (2, 3)), Word('x3', (2, 3))]
			>>> for x in missing: print(x)
			x2
			x3
		
		:rtype: a list of :class:`Words <thompson.word.Word>`."""
		words = set(self.simple_words_above())
		#At the end, should contract to [x1, x2, ..., x_r]
		expected = Generators.standard_basis(self.signature)
		missing = [x for x in expected if x not in words]
		return missing
	
	def generates_algebra(self):
		"""Returns True if this set generates all of :math:`V_{n,r}`. Otherwise returns False.
		
			>>> from random import randint
			>>> arity = randint(2, 5); alphabet_size = randint(2, 10)
			>>> Y = Generators.standard_basis((arity, alphabet_size))
			>>> Y.generates_algebra()
			True
			>>> random_basis().generates_algebra()
			True
			
		"""
		return len(self.test_generates_algebra()) == 0
		
	def test_above(self, word, return_index=False):
		"""Searches for a pair (gen, tail) where gen belongs to the current basis and gen + tail = word. If no such pair exists, returns None.
		
		If *return_index* is False, returns the pair *(gen, tail)*. Otherwise, returns the pair *(i, tail)* where *i* is the index of *gen* in the current basis.
		
			>>> basis = Generators.standard_basis((2, 1)).expand(0).expand(0).expand(0)
			>>> basis
			Generators((2, 1), ['x1 a1 a1 a1', 'x1 a1 a1 a2', 'x1 a1 a2', 'x1 a2'])
			>>> gen, tail = basis.test_above(Word('x1 a2 a2 a1', (2, 1)))
			>>> print(gen, '|', format(tail))
			x1 a2 | a2 a1
			>>> i, tail = basis.test_above(Word('x1 a2 a2 a1', (2, 1)), return_index=True)
			>>> print(basis[i])
			x1 a2
		"""
		if not return_index and self.cache is not None:
			return self._test_above_cached(word)
		for i, gen in enumerate(self):
			result = gen.test_above(word)
			if result is not None:
				if return_index:
					return i, result
				return gen, result
		return None
	
	def _test_above_cached(self, word):
		"""A quicker version of :meth:`test_above` which assumes that the:
			- a cache set has been maintained for this generating set, and
			- each generator is simple.
		
		.. todo::
			
			The bulk of the speedup seemed to be properly using all the information available. Do we still need this cache? Given the fact that the QNB is used virtually all the time it might be a good idea.
		"""
		if not word.is_simple():
			return None
		for i in range(len(word)+1):
			if word[:i] in self.cache:
				head, tail = word.split(i)
				return Word(head, word.signature, preprocess=False), tail
		return None
	
	def is_above(self, word):
		"""Returns True if any generator :meth:`~thompson.word.Word.is_above` the given *word*.
		
			>>> g = Generators((2, 2), ['x1 a1', 'x1 a2', 'x2'])
			>>> g.is_above(Word('x1 a1 a1 a2', (2, 2)))
			True
			>>> g.is_above(Word('x1', (2, 2)))
			False
			>>> g.is_above(Word('x2', (2, 2)))
			True
			>>> g.is_above(Word('x1 a1 x1 a2 L', (2, 2)))
			False
		"""
		return any(gen.is_above(word) for gen in self)
	
	def is_basis(self):
		"""Returns True if this is a :meth:`free generating set <is_free>` which :meth:`generates all of <is_free>` :math:`V_{n,r}`. Otherwise returns False.
		
			>>> g = Generators((2, 2), ["x1 a1", "x1 a2"])
			>>> g.is_free()
			True
			>>> g.generates_algebra()
			False
			>>> g.is_basis()
			False
			>>> g.append('x2')
			>>> g.is_basis()
			True
			>>> random_basis().is_basis()
			True
		"""
		return self.is_free() and self.generates_algebra()
	
	#Information about generating sets
	def descendants_above(self, floor):
		r"""Suppose we have a basis *floor* below the current basis *ceiling*. There are a finite number of elements below *ceiling* which are not below *floor*. This method enumerates them. In symbols, we are enumerating the set :math:`X\langle A\rangle \setminus F\langle A\rangle`, where :math:`X` is the current basis and :math:`F` is the *floor*.
		
			>>> X = example_4_25.quasinormal_basis
			>>> Y = X.minimal_expansion_for(example_4_25)
			>>> Z = example_4_25.image_of_set(Y)
			>>> terminal = X.descendants_above(Y)
			>>> initial  = X.descendants_above(Z)
			>>> print(initial, terminal, sep='\n')
			[x1 a1, x1 a1 a1]
			[x1 a2 a2, x1 a2 a2 a1]
		"""
		if not all(self.is_above(f) for f in floor):
			raise ValueError('Floor is not below the current generating set.')
		
		output = copy(self)
		i = len(output) - 1
		while i >= 0:
			if output[i] in floor:
				del output[i]
			i -= 1
		
		i = 0
		while i < len(output):
			word = output[i]
			children = [child for child in word.expand() if child not in floor]
			output.extend(children)
			i += 1
		return output
	
	#Creating generating sets
	@classmethod
	def standard_basis(cls, signature):
		r"""Creates the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` for :math:`V_{n,r}`, where :math:`n` is the arity and :math:`r` is the *alphabet_size* of the *signature*. 
		
			>>> Generators.standard_basis((2, 4))
			Generators((2, 4), ['x1', 'x2', 'x3', 'x4'])
		"""
		if isinstance(signature, tuple) and not isinstance(signature, Signature):
			signature = Signature(*signature)
		
		output = Generators(signature)
		for i in range(1, signature.alphabet_size + 1):
			output.append(Word((i,), signature, preprocess=False))
		return output
	
	def minimal_expansion_for(self, *automorphisms):
		r"""Suppose we are given a finite sequence of *automorphisms* of :math:`V_{n,r}` and that the current generaeting set is a basis :math:`X` for :math:`V_{n,r}`. This methods returns an expansion :math:`Y` of :math:`X` such that each automorphism maps :math:`Y` into :math:`X\langle A\rangle`.
		
			>>> std_basis = Generators.standard_basis((2, 1))
			>>> reduced_domain = Generators((2, 1), ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2", "x a2 a1", "x a2 a2"])
			>>> std_basis.minimal_expansion_for(cyclic_order_six) == reduced_domain
			True
			>>> basis = example_5_3.quasinormal_basis
			>>> print(basis.minimal_expansion_for(example_5_3))
			[x1 a1 a1 a1 a1, x1 a1 a1 a1 a2, x1 a1 a1 a2, x1 a1 a2 a1, x1 a1 a2 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2]
		
		:raises ValueError: if no automorphisms are passed to this method.
		:raises ValueError: if the automorphisms or basis do not all belong to the same group :math:`G_{n,r}`.
		:raises ValueError: if the given basis does not generate :math:`V_{n,r}`.
		
		.. seealso:: Lemma 4.2 of the paper proves that this expansion exists, is of minimal size, and is unique with this that property.
		"""
		#0. Pre-flight checks.
		if len(automorphisms) == 0:
			raise ValueError('Must provide at least one automorphism to minimal_expansion_for().')
		
		basis = copy(self)
		if not basis.is_free():
			raise ValueError('The generating set does not freely generate V_{n,r}.')
		
		for aut in automorphisms:
			if aut.signature != self.signature:
				raise ValueError('At least one automorphism belongs to a different G_{n,r} than the basis.')
		
		#1. Expand until all images belong to X<A>.
		i = 0
		while i < len(basis):
			b = basis[i]
			if all(self.is_above(aut.image(b)) for aut in automorphisms):
				i += 1
			else:
				basis.expand(i)
		return basis
	
	#Modifiers
	def expand(self, index):
		r"""Replaces the word :math:`w` at index *index* in the current generating set with its children, :math:`w\alpha1, \dotsc, w\alpha_n`, where :math:`n` is the arity of the generating set.
		
			>>> g = Generators.standard_basis((3, 1)); g
			Generators((3, 1), ['x1'])
			>>> g.expand(0)
			Generators((3, 1), ['x1 a1', 'x1 a2', 'x1 a3'])
			>>> g.expand(1)
			Generators((3, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2', 'x1 a2 a3', 'x1 a3'])
			
			:raises IndexError: if there is no generator at index *index*.
			:returns: the current generating set, after modification. 
		"""
		self[index: index+1] = [self[index].alpha(i) for i in range(1, self.signature.arity+1)]
		return self #allows chaining
	
	def _expand_with_cache(self, index):
		r"""If a cache is being maintained, this function expands() and maintains the cache."""
		assert self.cache is not None
		self.cache.remove(self[index])
		self.expand(index)
		self.cache.update(self[index : index + self.signature.arity])
		return self
	
	@classmethod
	def sort_mapping_pair(cls, domain, range):
		"""Makes copies of the given lists of words *domain* and *range*. The copy of *domain* is sorted according to the :meth:`order <thompson.word.Word.__lt__>` defined on words. The same re-ordering is applied to the *range*, so that the mapping from *domain* to *range* is preserved.
		
		:rtype: A pair of class:`Generator` instances.
		"""
		d, r = zip(*sorted(zip(domain, range)))
		return cls(domain.signature, d), cls(range.signature, r)
	
	def expand_to_size(self, size):
		"""Expands the current generating set until it has the given *size*. The expansions begin from the end of the generating set and work leftwards, wrapping around if we reach the start. (This is to try and avoid creating long words where possible.)
		
			>>> basis = Generators.standard_basis((3, 1)); print(basis)
			[x1]
			>>> basis.expand_to_size(11); print(basis)
			[x1 a1 a1, x1 a1 a2, x1 a1 a3, x1 a2 a1, x1 a2 a2, x1 a2 a3, x1 a3 a1, x1 a3 a2, x1 a3 a3 a1, x1 a3 a3 a2, x1 a3 a3 a3]
			>>> basis = Generators.standard_basis((2, 1)); print(basis)
			[x1]
			>>> basis.expand_to_size(2); print(basis)
			[x1 a1, x1 a2]
			>>> basis = Generators.standard_basis((2, 3)).expand(2).expand(0); print(basis)
			[x1 a1, x1 a2, x2, x3 a1, x3 a2]
			>>> basis.expand_to_size(12); print(basis)
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2, x2 a1, x2 a2, x3 a1 a1, x3 a1 a2, x3 a2 a1 a1, x3 a2 a1 a2, x3 a2 a2 a1, x3 a2 a2 a2]
		
		:returns: the current generating set
		
		:raises ValueError: if an expansion to the given size is not possible.
		"""
		modulus = self.signature.arity - 1
		if (size % modulus != len(self) % modulus) or size < len(self):
			raise ValueError("Cannot expand from length {} to length {} in steps of size {}.".format(
			  len(self), size, modulus))
			  
		num_expansions = (size - len(self))// modulus
		i = -1
		for _ in range(num_expansions):
			if i < 0:
				i = len(self) - 1
			self.expand(i)
			i -= 1
		assert len(self) == size

def rewrite_set(set, basis, new_basis):
	"""Maps *set* into the algebra generated by *new_basis* according to the bijection between *basis* and *new_basis*. See also :func:`rewrite_word`."""
	set.signature = new_basis.signature
	for i, word in enumerate(set):
		set[i] = rewrite_word(word, basis, new_basis)

def rewrite_word(word, basis, new_basis):
	r"""Suppose that we have a *word* which is below a given *basis*. Suppose we have a bijection between *basis* and some image *new_basis*. Then we can rewrite *word* as as descendant of *new_basis* in a manner which is compatible with this bijection.
	
	:raises ValueError: if *basis* is not above *word*.
	:raises ValueError: if *basis* and *new_basis* have different sizes.
	"""
	if len(basis) != len(new_basis):
		raise ValueError('Basis sizes {} and {} do not match.'.format(
		  len(basis), len(new_basis)))
	
	result = basis.test_above(word, return_index=True)
	if result is None:
		raise ValueError('The word {} is not below the basis {}'.format(
		  word, basis))
	index, tail = result
	return new_basis[index].extend(tail)
