from copy import copy

from .word import Word

class Generators(list):
	"""A list of :class:`~thompson.word.Words`, together with methods which treat this as a Basis for :math:`V_{n,r}`. The list is assumed to not contain any duplicates."""
	
	def __init__(self, *words):
		"""Creates a generator object from its arguments. Each argument *word* in *words* should be either a :class:`~thompson.word.Word` or a :meth:`string describing one <thompson.word.Words.from_string>`.
		
		:raises TypeError: if one of the *words* is neither a string nor a Word.
		:raises ValueError: if *words* contains duplicated words."""
		
		seen = set()
		parsed = list(words)
		
		for i, word in enumerate(parsed):
			if isinstance(word, str):
				word = parsed[i] = Word.from_string(word)
			elif not isinstance(word, Word):
				raise TypeError("Argument {:r} at position {} is neither a string nor a Word.".format(word, i))
			representation = str(word)
			if representation not in seen:
				word.standardise()
				seen.add(representation)
			else:
				raise ValueError("The word {} at position {} is a duplicate.".format(word, i))
		
		super().__init__(self)
		self.extend(parsed)
	
	def is_free(self):
		"""Returns True if this set of generators is free; otherwise False."""
		return self.test_free() == (-1, -1)
	
	def test_free(self):
		"""
		Tests to see if this set of generators is free. If the test fails, returns the first pair of indices (i, j) found for which one of self[i] and self[j] is an initial segment of the other. If the test passes, returns (-1, -1) 
		
		Lemma 3.16.3: :math:`V` is a set of *free* generators for the subalgebra that it generates iff no element of :math:`V` is a proper initial segment of another.
		"""
		for i in range(len(self)):
			for j in range(i+1, len(self)):
				#Assumption: self contains no duplicates 
				if Word.initial_segment(self[i], self[j]):
					return i, j
		return -1, -1
	
	def test_basis(self, arity, alphabet_size):
		"""Returns True if this set is an expansion of :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r`\}`, where :math:`n` is *arity* and :math:`r` is the *alphabet_size*. Otherwise, returns false."""
		words = sorted(copy(self)) #shallow copy: read only
		words.sort()
		#This part is very similar to :meth:`thompson.word.Word._contract`.
		i = 0
		while i < len(words) - arity:
			prefix = Word.are_contractible(words[i : i + arity])
			if prefix is not None:
				words[i : i + arity] = [prefix]
			else:
				i += 1
		
		#At the end, should contract to [x1, x2, ..., x_r]
		return all(word == [i + 1] for i, word in enumerate(words))
