
	def __init__(self, words, arity=None, alphabet_size=None):
		"""A set of generators for :math:`V_{n, r}`.
		
		:param words: A list of words. Each entry should be either a :class:`~thompson.word.Word` or a :meth:`string describing one <thompson.word.Words.from_string>`.
		:param arity: If given
		
		:raises TypeError: if one of the *words* is neither a string nor a Word.
		:raises ValueError: if one of the *words* has a different arity to *arity*.
		:raises ValueError: if one of the *words* has a larger alphabet_size than *alphabet_size*.
		:raises ValueError: if *words* contains duplicated words.
		
			>>> Generators([1234, 'x a2'], 2, 1)
				...
			TypeError(Word 1234 at position 0 is neither a string nor a Word.)
			>>> Generators(['x a1', 'x a1'], 2, 1)
			
		"""
	
	def is_free(self):
		"""Returns True if this is a free generating set; otherwise False."""
		return self.test_free() == (-1, -1)
	
	def test_free(self):
		"""Tests to see if this is a free generating set. If the test fails, returns the first pair of indices (i, j) found for which one of self[i] and self[j] is an initial segment of the other. If the test passes, returns (-1, -1). 
		
		Lemma 3.16.3: :math:`V` is a set of *free* generators for the subalgebra that it generates iff no element of :math:`V` is a proper initial segment of another.
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
		#This part is very similar to :meth:`~thompson.word.Word._contract`.
		i = 0
		while i < len(words) - self.arity:
			prefix = Word.are_contractible(words[i : i + self.arity])
			if prefix is not None:
				words[i : i + self.arity] = [prefix]
			else:
				i += 1
		
		#At the end, should contract to [x1, x2, ..., x_r]
		return all(word == [i + 1] for i, word in enumerate(words))
	
	@classmethod
	def standard_basis(cls, arity):
		"""Returns the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r`\}` of :math:`V_{n,r}`, where :math`r` is the *arity*."""
		generators = (Word().letter(i) for i in range(1, arity+1))
		return cls(*generators)
	