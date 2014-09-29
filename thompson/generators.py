"""
.. testsetup:: 
	
	from thompson.generators import *
"""
from .word import Word, are_contractible

__all__ = ["Generators"]

class Generators(list):
	r"""An ordered subset of :math:`V_{n, r}`, together with methods which treat such sets as generating sets and bases. Internally, this is a subclass of :class:`list <py3:list>`, so it is possible to reorder the generating set.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	"""
	
	def __init__(self, arity, alphabet_size):
		"""When creating a generating set, you must specify the algebra :math:`V_{n, r}` it is a subset of."""
		self.arity = arity
		self.alphabet_size = alphabet_size
	
	def append(self, word):
		"""Adds *word* to this generating set. If *word* is a string, a :class:`Word` object is created and assigned the same arity and alphabet_size as the generating set.
		
		:raises TypeError: if *word* is neither a string nor a Word.
		:raises ValueError: if one of the *words* has a different arity to the generating set.
		:raises ValueError: if one of the *words* has a larger alphabet_size the generating set.
		:raises ValueError: if *word* is already contained in this generating set.
		"""
		if isinstance(word, str):
			word = Word(word, self.arity, self.alphabet_size)
		elif not isinstance(word, Word):
			raise TypeError("{:r} is neither a string nor a Word.".format(word))
		
		if word.arity != self.arity:
			raise ValueError("Can't add {} with arity {} to generating set with arity {}."
			  .format(word, word.arity, self.arity))
		
		if word.alphabet_size > self.alphabet_size:
			raise ValueError("Can't add {} with alphabet size {} to generating set with alphabet size {}.".format(
			  word, word.alphabet_size, self.alphabet_size))
		#else: modify the alphabet size?
		
		if word not in self:
			super().append(word)
		else:
			raise ValueError("The word {} is already in this set.".format(word))
		
	def is_free(self):
		"""Returns True if this is a free generating set; otherwise False.
		
			>>> g = Generators(2, 3)
			>>> g.append('x1'); g.append('x2');
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
		for i in range(len(self)):
			for j in range(i+1, len(self)):
				#Assumption: self contains no duplicates 
				if Word.initial_segment(self[i], self[j]):
					return i, j
		return -1, -1
	
	def is_basis(self):
		"""Returns True if this set is an expansion of :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r`\}`."""
		words = sorted(self) #Creates a shallow copy and then sorts it.
		#This part is very similar to :meth:`~thompson.word._contract`.
		i = 0
		while i < len(words) - self.arity:
			prefix = are_contractible(words[i : i + self.arity])
			if prefix is not None:
				words[i : i + self.arity] = [prefix]
			else:
				i += 1
		
		#At the end, should contract to [x1, x2, ..., x_r]
		return all(word == [i + 1] for i, word in enumerate(words))
	
	@classmethod
	def standard_basis(cls, arity, alphabet_size):
		"""Creates the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r`\}` for :math:`V_{n,r}`, where :math:`n` is the arity and :math`r` is the *alphabet_size*."""
		generators = Generators(arity, alphabet_size)
		for i in range(alphabet_size):
			generators.append(Word([i+1], arity, alphabet_size))
		return generators

