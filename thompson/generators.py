"""
.. testsetup:: 
	
	from thompson.generators import *
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
		:raises ValueError: if *w* has a different arity to the generating set.
		:raises ValueError: if *w* has a larger alphabet_size the generating set.
		:raises ValueError: if *w* is already contained in this generating set.
		"""
		if isinstance(w, str):
			w = Word(w, self.arity, self.alphabet_size)
		elif not isinstance(w, Word):
			raise TypeError("{:r} is neither a string nor a Word.".format(w))
		
		if w.arity != self.arity:
			raise ValueError("Can't add {} with arity {} to generating set with arity {}."
			  .format(w, w.arity, self.arity))
		
		if w.alphabet_size > self.alphabet_size:
			raise ValueError("Can't add {} with alphabet size {} to generating set with alphabet size {}.".format(
			  w, w.alphabet_size, self.alphabet_size))
		#else: modify the alphabet size?
		
		if w not in self:
			super().append(w)
		else:
			raise ValueError("The word {} is already in this set.".format(w))
		
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
	
	def test_free(self):
		"""Tests to see if this is a free generating set. If the test fails, returns the first pair of indices (i, j) found for which one of self[i] and self[j] is an initial segment of the other. If the test passes, returns (-1, -1). See lemma 3.16.3.
		"""
		#TODO example
		for i in range(len(self)):
			for j in range(i+1, len(self)):
				#Assumption: self contains no duplicates 
				if word.initial_segment(self[i], self[j]):
					return i, j
		return -1, -1
	
	def test_generates_algebra(self):
		"""Tests to see if this set generates all of :math:`V_{n,r}`. The generating set is sorted, and then contracted as much as possible. Then we check to see which elements of the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` are not included in the contraction.
		
		The test fails if there is at least one such element; in this case, the list of all such elements is returned. (The elements of this list are lists of integers, not fully-fledged :class:`Words <Word>`.) Otherwise the test passes, and an empty list is returned.
		"""
		#Todo doctest.
		words = sorted(self) #Creates a shallow copy and then sorts it.
		#This part is very similar to :meth:`~thompson.word._contract`.
		i = 0
		while i < len(words) - self.arity:
			prefix = word.are_contractible(words[i : i + self.arity])
			if prefix is not None:
				words[i : i + self.arity] = [prefix]
			else:
				i += 1
		
		#At the end, should contract to (something which includes) [x1, x2, ..., x_r]
		missing = []
		for i in range(self.alphabet_size):
			if [i] not in words:
				missing.append([i])
		return missing
	
	@classmethod
	def standard_basis(cls, arity, alphabet_size):
		"""Creates the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r`\}` for :math:`V_{n,r}`, where :math:`n` is the arity and :math`r` is the *alphabet_size*."""
		#Todo example.
		generators = Generators(arity, alphabet_size)
		for i in range(alphabet_size):
			generators.append(Word([i+1], arity, alphabet_size))
		return generators
	
	#todo __str__

