"""
.. testsetup:: 
	
	from thompson.generators import *
	from thompson import word
"""
from . import word
from .word import Word

__all__ = ["Generators"]

class Generators(list):
	r"""An ordered subset of :math:`V_{n, r}`, together with methods which treat such sets as generating sets and bases. Internally, this is a subclass of :class:`list <py3:list>`, so it is possible to reorder the generating set.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	"""
	
	def __init__(self, arity, alphabet_size, generators=None):
		"""When creating a generating set, you must specify the algebra :math:`V_{n, r}` it is a subset of. A list of generators can also be provided to save calls to :meth:`append`."""
		self.arity = arity
		self.alphabet_size = alphabet_size
		if generators is not None:
			for g in generators:
				self.append(g)
	
	def append(self, w):
		"""Adds *w* to this generating set. If *w* is a string, a :class:`Word` object is created and assigned the same arity and alphabet_size as the generating set.
		
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
	
	def test_free(self):
		"""Tests to see if this is a free generating set. If the test fails, returns the first pair of indices (i, j) found for which one of ``self[i]`` and ``self[j]`` is an initial segment of the other. If the test passes, returns (-1, -1). See lemma 3.16.3.
		
			>>> Generators(2, 3, ["x1 a1", "x1 a2 a1", "x1 a2 a1 a1", "x1 a2 a2"]).test_free()
			(1, 2)
		
		"""
		for i in range(len(self)):
			for j in range(i+1, len(self)):
				#Assumption: self contains no duplicates 
				if word.initial_segment(self[i], self[j]):
					return i, j
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
		"""Tests to see if this set generates all of :math:`V_{n,r}`. The generating set is sorted, and then contracted as much as possible. Then we check to see which elements of the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` are not included in the contraction.
		
		The test fails if there is at least one such element; in this case, the list of all such elements is returned. (The elements of this list are lists of integers, not fully-fledged :class:`Words <Word>`.) Otherwise the test passes, and an empty list is returned.
		
			>>> #A basis for V_2,1
			>>> Y = Generators(2, 1, ["x a1", "x a2 a1", "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
			>>> Y.test_generates_algebra()
			[]
			>>> #The same words do not generate V_2,3
			>>> Y = Generators(2, 3, Y)
			>>> missing = Y.test_generates_algebra(); missing
			[[2], [3]]
			>>> print(", ".join(word.format(x) for x in missing))
			x2, x3
		"""
		words = sorted(self) #Creates a shallow copy and then sorts it.
		#This part is very similar to :meth:`~thompson.word._contract`.
		i = 0
		while i <= len(words) - (self.arity - 1):
			prefix = word.are_contractible(words[i : i + self.arity])
			if prefix is not None:
				words[i : i + self.arity] = [prefix]
				i -= (self.arity - 1) 
				i = max(i, 0)
			else:
				i += 1
		
		#At the end, should contract to (something which includes) [x1, x2, ..., x_r]
		missing = []
		for i in range(1, self.alphabet_size + 1):
			if (i,) not in words:
				missing.append([i])
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
	
	@classmethod
	def standard_basis(cls, arity, alphabet_size):
		"""Creates the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r`\}` for :math:`V_{n,r}`, where :math:`n` is the arity and :math`r` is the *alphabet_size*. 
		
			>>> Generators.standard_basis(2, 4)
			Generators(2, 4, ['x1', 'x2', 'x3', 'x4'])
		"""
		output = Generators(arity, alphabet_size)
		for i in range(1, alphabet_size + 1):
			output.append(Word([i], arity, alphabet_size))
		return output
	
	def expand(self, index):
		r"""Replaces the word :math`w` at index *index* with its children, :math:`w\alpha1, \dotsc, w\alpha_n`, where :math`n` is the arity of the generating set.
		
			>>> g = Generators.standard_basis(3, 1); g
			Generators(3, 1, ['x1'])
			>>> g.expand(0); g
			Generators(3, 1, ['x1 a1', 'x1 a2', 'x1 a3'])
			>>> g.expand(1); g
			Generators(3, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2', 'x1 a2 a3', 'x1 a3'])
			
			:raises IndexError: if there is no generator at index *index*.
		"""
		self[index: index+1] = [self[index].alpha(i) for i in range(1, self.arity+1)]

