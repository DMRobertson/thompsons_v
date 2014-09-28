from .word import Word

class Generators:
	r"""An ordered subset of :math:`V_{n, r}`, together with methods which treat such sets as generating sets and bases.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	"""
	
	def __init__(self, arity, alphabet_size):
		"""When creating a generating set, you must specify the algebra :math:`V_{n, r}` it belongs to."""
		self.arity = arity
		self.alphabet_size = alphabet_size
		self._words = []
	
	def append(self, word):
		"""Adds *word* to this generating set. If *word* is a string, it is passed to :meth:`~thompson.word.Word.from_string`to form a :class:`Word`.
		
		:raises TypeError: if *word* is neither a string nor a Word.
		:raises ValueError: if one of the *words* has a different arity to the generating set.
		:raises ValueError: if one of the *words* has a larger alphabet_size the generating set.
		:raises ValueError: if *word* is already contained in this generating set.
		"""
		if isinstance(word, str):
			word = Word.from_string(word)
		elif not isinstance(word, Word):
			raise TypeError("{:r} is neither a string nor a Word.".format(word))
		
		if word.arity != self.arity:
			raise ValueError("Can't add {} with arity {} to generating set with arity {}."
			  .format(word, word.arity, self.arity))
		
		if word.alphabet_size > alphabet_size:
			raise ValueError("Can't add {} with alphabet size {} to generating set with alphabet size {}.".
			  .format(word, word.alphabet_size, alphabet_size))
		
		if word not in self._words:
			self._words.append(word)
		else:
			raise ValueError("The word {} is already in this set.".format(word))
		
		