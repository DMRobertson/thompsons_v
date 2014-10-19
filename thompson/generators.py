r"""Let :math:`S` be a set of words in :math:`V_{n,r}`. We can obtain new words by applying the operations :math:`\alpha_i` and :math:`\lambda` to elements of `S`. Next, we can apply the operations to the new words we've received. This process continues indefinitely, producing a set of words :math:`T` whose elements all come from the original set :math:`S`. We call :math:`T` the *subalgebra generated by* :math:`S`.

Of particular interest are *bases*, which we can think of as generating sets which do not contain any redundant elements. For example, it would be pointless to include each of :math:`x`, :math:`x \alpha_1` and :math:`x \alpha_2` in `S`, as the first can be obtained from the last two and vice versa.

.. testsetup:: 
	
	from thompson.generators import *
	from thompson.word import Word, format

The Generators class
--------------------
"""
from copy import copy

from . import word
from .word import Word

__all__ = ["Generators"]

class Generators(list):
	r"""An ordered subset of :math:`V_{n, r}`, together with methods which treat such sets as generating sets and bases. Internally, this is a subclass of :class:`list <py3:list>`, so it is possible to reorder and otherwise modify collections of Generators.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	"""
	
	def __init__(self, arity, alphabet_size, generators=None):
		"""When creating a generating set, you must specify the algebra :math:`V_{n, r}` it is a subset of. A list of generators can optionally be provided. Each generator is passed to :meth:`append`."""
		self.arity = arity
		self.alphabet_size = alphabet_size
		if generators is not None:
			for g in generators:
				self.append(g)
	
	def append(self, w):
		"""Adds *w* to this generating set. If *w* is a string, a :class:`Word` object is created and assigned the same *arity* and *alphabet_size* as the generating set.
		
		:raises TypeError: if *w* is neither a string nor a Word.
		:raises IndexError: if *w* has a different arity to the generating set.
		:raises IndexError: if *w* has a larger alphabet_size the generating set.
		:raises ValueError: if *w* is already contained in this generating set.
		"""
		if isinstance(w, str):
			w = Word(w, self.arity, self.alphabet_size)
		elif not isinstance(w, Word):
			raise TypeError("{:r} is neither a string nor a Word.".format(w))
		
		if w.arity != self.arity:
			raise IndexError("Can't add {} with arity {} to generating set with arity {}."
			  .format(w, w.arity, self.arity))
		
		if w.alphabet_size > self.alphabet_size:
			raise IndexError("Can't add {} with alphabet size {} to generating set with alphabet size {}.".format(
			  w, w.alphabet_size, self.alphabet_size))
		#else: modify the alphabet size?
		
		if w not in self:
			super().append(w)
		else:
			raise ValueError("The word {} is already in this set.".format(w))
	
	def __str__(self):
		return "[" + ", ".join(format(w) for w in self) + "]"
	
	def __repr__(self):
		return "Generators({}, {}, [{}])".format(
		  self.arity, self.alphabet_size, ", ".join(repr(format(w)) for w in self))
	
	#Tests on generating sets
	def test_free(self):
		r"""Tests to see if the current generating set ``X`` is a free generating set. If the test fails, returns the first pair of indices :math:`(i, j)` found for which one of ``X[i]`` and ``X[j]`` is an initial segment of the other. If the test passes, returns ``(-1, -1)``.
		
			>>> g = Generators(2, 3, ["x1 a1", "x1 a2 a1", "x1 a2 a1 a1", "x1 a2 a2"])
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
		
			>>> g = Generators(2, 3, ['x1', 'x2']);
			>>> g.is_free()
			True
			>>> g.append('x3'); g.is_free()
			True
			>>> g.append('x2 a1'); g.is_free()
			False
		"""
		return self.test_free() == (-1, -1)
	
	def test_generates_algebra(self):
		r"""Tests to see if this set generates all of :math:`V_{n,r}`. The generating set is sorted, and then contracted as much as possible. Then we check to see which elements of the :meth:`standard basis <standard_basis>` :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` are not included in the contraction.
		
		The test fails if there is at least one such element; in this case, the list of all such elements is returned. Otherwise the test passes, and an empty list is returned.
		
			>>> #A basis for V_2,1
			>>> Y = Generators(2, 1, ["x a1", "x a2 a1", "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
			>>> Y.test_generates_algebra()
			[]
			>>> #The same words do not generate V_2,3
			>>> Y = Generators(2, 3, Y)
			>>> missing = Y.test_generates_algebra(); missing
			[Word('x2', 2, 3), Word('x3', 2, 3)]
			>>> for x in missing: print(x)
			x2
			x3
		
		:rtype: a list of :class:`Words <thompson.word.Word>`."""
		
		words = sorted(self) #Creates a shallow copy and then sorts it.
		#Todo maybe this should be part of a generators.contract method
		i = 0
		while i <= len(words) - self.arity:
			prefix = word.are_contractible(words[i : i + self.arity])
			if prefix:
				words[i : i + self.arity] = [prefix]
				i -= (self.arity - 1) 
				i = max(i, 0)
			else:
				i += 1
		
		#At the end, should contract to (something which includes) [x1, x2, ..., x_r]
		expected = Generators.standard_basis(self.arity, self.alphabet_size)
		missing = [x for x in expected if x not in words]
		return missing
	
	def generates_algebra(self):
		"""Returns True if this set generates all of :math:`V_{n,r}`. Otherwise returns False.
		
			>>> from random import randint
			>>> arity = randint(2, 5); alphabet_size = randint(2, 10)
			>>> Y = Generators.standard_basis(arity, alphabet_size)
			>>> Y.generates_algebra()
			True
		"""
		return len(self.test_generates_algebra()) == 0
		
	def test_above(self, word):
		"""Returns a pair (gen, tail) where gen belongs to the current basis and gen + tail = word. If no such pair exists, returns None.
		
			>>> basis = Generators.standard_basis(2, 1).expand(0).expand(0).expand(0)
			>>> basis
			Generators(2, 1, ['x1 a1 a1 a1', 'x1 a1 a1 a2', 'x1 a1 a2', 'x1 a2'])
			>>> gen, tail = basis.test_above(Word('x1 a2 a2 a1', 2, 1))
			>>> print(gen, '|', format(tail))
			x1 a2 | a2 a1
		"""
		for gen in self:
			result = gen.test_above(word)
			if result is not None:
				return gen, result
		return None
	
	def is_above(self, word):
		"""Returns True if any generator :meth:`~thompson.word.Word.is_above` the given *word*.
		
			>>> g = Generators(2, 2, ['x1 a1', 'x1 a2', 'x2'])
			>>> g.is_above(Word('x1 a1 a1 a2', 2, 2))
			True
			>>> g.is_above(Word('x1', 2, 2))
			False
			>>> g.is_above(Word('x2', 2, 2))
			True
			>>> g.is_above(Word('x1 a1 x1 a2 L', 2, 2))
			False
		"""
		return any(gen.is_above(word) for gen in self)
	
	def is_basis(self):
		"""Returns True if this is a :meth:`free generating set <is_free>` which :meth:`generates all of <is_free>` :math:`V_{n,r}`. Otherwise returns False.
		
			>>> g = Generators(2, 2, ["x1 a1", "x1 a2"])
			>>> g.is_free()
			True
			>>> g.generates_algebra()
			False
			>>> g.is_basis()
			False
			>>> g.append('x2')
			>>> g.is_basis()
			True
		"""
		return self.is_free() and self.generates_algebra()
	
	#Creating generating sets
	@classmethod
	def standard_basis(cls, arity, alphabet_size):
		r"""Creates the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` for :math:`V_{n,r}`, where :math:`n` is the arity and :math:`r` is the *alphabet_size*. 
		
			>>> Generators.standard_basis(2, 4)
			Generators(2, 4, ['x1', 'x2', 'x3', 'x4'])
		"""
		output = Generators(arity, alphabet_size)
		for i in range(1, alphabet_size + 1):
			output.append(Word((i,), arity, alphabet_size, preprocess=False))
		return output
	
	def minimal_expansion_for(self, *automorphisms):
		r"""Suppose we are given a finite sequence of *automorphisms* of :math:`V_{n,r}` and that the current generaeting set is a basis :math:`X` for :math:`V_{n,r}`. This methods returns an expansion :math:`Y` of :math:`X` such that each automorphism maps :math:`Y` into :math:`X\langle A\rangle`.

		
			>>> from thompson.examples import *
			>>> std_basis = Generators.standard_basis(2, 1)
			>>> std_basis.minimal_expansion_for(example_4_5) == example_4_5.domain
			True
			>>> std_basis.minimal_expansion_for(example_4_11) == example_4_11.domain
			True
			>>> qnf_basis = example_4_12.quasinormal_basis()
			>>> qnf_basis.minimal_expansion_for(example_4_12) == qnf_basis
			True
			>>> std_basis.minimal_expansion_for(example_4_25) == example_4_25.domain
			True
			>>> reduced_domain = Generators(2, 1, ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2", "x a2 a1", "x a2 a2"])
			>>> std_basis.minimal_expansion_for(cyclic_order_six) == reduced_domain
			True
			>>> basis = example_5_3.quasinormal_basis()
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
			if not (aut.arity == self.arity and aut.alphabet_size == self.alphabet_size):
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
		
			>>> g = Generators.standard_basis(3, 1); g
			Generators(3, 1, ['x1'])
			>>> g.expand(0)
			Generators(3, 1, ['x1 a1', 'x1 a2', 'x1 a3'])
			>>> g.expand(1)
			Generators(3, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2', 'x1 a2 a3', 'x1 a3'])
			
			:raises IndexError: if there is no generator at index *index*.
			:returns: the current generating set, after modification. 
		"""
		self[index: index+1] = [self[index].alpha(i) for i in range(1, self.arity+1)]
		return self #allows chaining

