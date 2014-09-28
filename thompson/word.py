r"""The algebra :math:`V_{n, r}` consists of words written in the alphabet :math:`\{x_1,\dotsc, x_r, \alpha_1, \dotsc, \alpha_n, \lambda\}`, as defined in remark 3.3. In the workings of this module, we represent a word as a list (or any iterable) of integers, where:

- :math:`x_i` is represented by :math:`i`,
- :math:`\alpha_i` is represented by :math:`-i`, and
- :math:`\lambda` is represented by :math:`0`.

.. testsetup:: 
	
	from thompson.word import *
"""
from itertools import chain

__all__ = ["format", "validate", "standardise", "are_contractible", "from_string", "Word"]

def _char(symbol):
	if symbol > 0:
		return 'x' + str(symbol)
	elif symbol < 0:
		return 'a' + str(-symbol)
	elif symbol == 0:
		return 'L'

def validate(letters, arity, alphabet_size=float('inf')):
	r"""Checks that the given list of *letters* makes sense as a word. If no errors are raised when running this function, then we know that *letters* represents a valid word. Specifically, we know that:
	
	- The first letter is an x.
	- Every :math:`\lambda`-contraction has *arity* arguments to its left.
	- No :math:`\alpha_i` appears with :math:`i >` *arity*.
	- No :math:`x_i` occurs with  :math:`i >` *alphabet_size*. This last check occurs only when *alphabet_size* is given.
	
	The *letters* may not describe a word in standard form---for that, use :func:`standardise`.
	
		>>> validate(from_string('a1 x a2 L'), arity=2)
		Traceback (most recent call last):
			...
		ValueError: The first letter (a1) should be an x.
		>>> validate(from_string('x1 a2 x12 a1 a2 L'), arity=2, alphabet_size=4)
		Traceback (most recent call last):
			...
		ValueError: Letter x12 at index 2 is not in the alphabet (maximum is x4).
		>>> validate(from_string('x1 a1 a2 a3 a4 a5 x2 a2 L'), arity=3)
		Traceback (most recent call last):
			...
		ValueError: Letter a4 at index 4 is not in the alphabet (maximum is a3).
		>>> validate(from_string('x1 a2 L'), arity=4)
		Traceback (most recent call last):
			...
		ValueError: Letter L at index 2 is invalid. Check that lambda has 4 arguments.
		>>> validate(from_string('x1 x2 L x3 L x4'), arity=2)
		Traceback (most recent call last):
			...
		ValueError: Word is invalid: valency is 2 (should be 1).
	
	:raises ValueError: if this word fails the valency test (prop 2.12).
	:raises ValueError: if this word contains an :math:`x_i` outside of the range 1 ... *alphabet_size*
	:raises ValueError: if this word contains an :math:`\alpha_i` outside of the range 1 ... *arity*
	"""
	symbol = letters[0]
	if symbol < 0:
		raise ValueError("The first letter ({}) should be an x.".format(
		  _char(symbol)))
	valency = 0
	for i, symbol in enumerate(letters):
		if symbol > alphabet_size:
			raise ValueError("Letter {} at index {} is not in the alphabet (maximum is x{}).".format(
			  _char(symbol), i, alphabet_size))
		if symbol < -arity:
			raise ValueError("Letter {} at index {} is not in the alphabet (maximum is a{}).".format(
			  _char(symbol), i, arity))
		valency += _valency_of(symbol, arity)
		if valency <= 0:
			raise ValueError("Letter {} at index {} is invalid. Check that lambda has {} arguments.".format(
			  _char(symbol), i, arity))
	if valency != 1:
		raise ValueError("Word is invalid: valency is {} (should be 1).".format(valency))

def _valency_of(symbol, arity):
	"""Returns the valency of a symbol as per Def 2.11."""
	if symbol > 0: #letter x_i
		return 1
	if symbol < 0: #descenders alpha
		return 0
	if symbol == 0:#contractor lambda
		return 1 - arity

def standardise(letters, arity):
	"""Accepts a word (as a list of letters) and reduces it to standard form. The result (a new list of letters) is then returned. See also remark 3.3.
	
	In the examples below, the :class:`Word` class standardises its input before creating a Word.
	
		>>> #Extraction
		>>> print(Word("x a1 a2 a1 x a2 a1 L a1 a2", 2, 1))
		x1 a1 a2 a1 a2
		>>> #Contraction
		>>> print(Word("x2 a2 a2 a1 x2 a2 a2 a2 L", 2, 2))
		x2 a2 a2
		>>> print(Word("x1 a1 a1 x1 a1 a2 x1 a1 a3 L", 3, 2))
		x1 a1
		>>> #Both extraction and contraction
		>>> print(Word("x a1 a2 a1 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L", 2, 1))
		x1 a1 a2 a1
		>>> #Something that doesn't simplify all the lambdas away
		>>> print(Word("x a1 a2 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L", 2, 1))
		x1 a1 a2 a1 x1 a1 a2 a1 a2 L
	
	:raises IndexError: if the list of letters :func:`is not valid <validate>`..
	"""
	#Convert the iterable to a list so that we can modify it in place.
	letters = list(letters)
	subwords = []
	i = 0
	while i < len(letters):
		symbol = letters[i]
		if symbol > 0: #an x -> new subword
			subwords.append([symbol])
			i += 1
		elif symbol < 0: #an alpha
			subwords[-1].append(symbol)
			i += 1
		else: #lambda. Run validate() to ensure it has enough arguments.
			#Peek ahead: is the lambda immediately followed by an alpha?
			try:
				next = letters[i+1]
			except IndexError:
				next = 0
			
			if next < 0:
				i = _extract(letters, subwords, arity, i, -next)
			else:
				i = _contract(letters, subwords, arity, i)
		#end lambda
		
	return letters

def _extract(letters, subwords, arity, lambda_position, index):
	r"""If we find a substring of the form :math`w_1 \dots w_\text{arity} \lambda \alpha_\text{index}`, replace it by :math:`w_\text{index}`.
	
	:return: Returns the index in *letters* to move to after extraction."""
	start = lambda_position - sum(len(list) for list in subwords[-arity:])
	extracted = subwords[-arity + index - 1]
	
	del subwords[-arity + 1:]
	subwords[-1] = extracted
	letters[start : lambda_position + 2] = extracted
	return start + len(extracted)

def are_contractible(words):
	r"""Let *words* be a list of words, either as lists of integers or as full :class:`Word` objects. This function tests to see if *words* is a list of the form :math:`(w\alpha_1, \dotsc, w\alpha_n)`, where ``n == len(words)``.
	
	:return: :math:`w` (as a tuple of integers) if the test passes; ``None`` if the test fails.
	
	**NB:** If there is a prefix :math:`w`, this function does **not** check to see if x :func:`is a valid word <validate>`.
	
		>>> are_contractible([Word("x a1 a2 a3 a1", 3, 1), Word("x a1 a2 a3 a2", 3, 1), Word("x a1 a2 a3 a3", 3, 1)])
		(1, -1, -2, -3)
	"""
	#TODO Slicing a Word gives a tuple. One the one hand, this is saying that "this won't neccesarily be a word", and Word objects are supposed to be fixed in standard form. That makes sense. But it would be convenient for this example to return a word. Then again, this function knows nothing of arity and alphabet size, so it can't generate a word.
	#Probably best to export word.__str__ to a module function and use this in the example.
	prefix = words[0][:-1]
	if len(prefix) == 0:
		return None
	expected_length = len(prefix) + 1
	
	for j, word in enumerate(words):
		if not (
		  len(word) == expected_length
		  and word[:len(prefix)] == prefix
		  and word[-1] == -j - 1): #alpha_{j+1}
			return None
	return prefix

def _contract(letters, subwords, arity, lambda_position):
	"""Attempts to contract the last *arity* words in *subwords*.
	
	Returns the index of the index in *letters* to move to after contraction.
	"""
	inspected = subwords[-arity:]
	start = lambda_position - sum(len(list) for list in inspected)
	prefix = are_contractible(inspected)
	
	if prefix is not None:
		#TODO: The lambda call is not of the form ua1, ua2, ..., uan, lambda, and the next symbol is NOT an alpha.
		#I think there's no way this can be reduced further.
		subwords[-arity] = prefix
		del subwords[-arity + 1:]
		letters[start : lambda_position + 1] = prefix
		return start + len(prefix)
	else:
		#We can contract the last *arity* words.
		subwords[-arity] = list(chain.from_iterable(inspected))
		del subwords[-arity + 1:]
		return lambda_position + 1 #plus one: can't do anything to this lambda, move on

def from_string(str):
	"""Converts a string representation of a word to the internal format (a list of integers). Anything which does not denote a basis letter (e.g. ``'x'``, ``'x2'``, ``'x45'``) a descendant operator (e.g.  ``'a1'``, ``'x3'``, ``'x27'``) or a lambda contraction (``'L'``) is ignored.
	"""
	output = str.lower().split()
	for i, string in enumerate(output):
		char = string[0]
		if char == 'x':
			value = 1 if len(string) == 1 else int(string[1:])
			output[i] = value
		elif char == 'a':
			value = int(string[1:])
			output[i] = -value
		elif char == ('l'):
			output[i] = 0
	return output

def format(letters):
	"""Turns a list of letters into a string. Can be thought of as an inverse to :func:`from_string`."""
	return " ".join(_char(i) for i in letters)

class Word(tuple):
	r"""Represents an element of the algebra :math:`V_{n, r}`. This consists of words written in the alphabet :math:`\{x_1,\dotsc, x_r, \alpha_1, \dotsc, \alpha_n, \lambda\}`---see Remark 3.3.
	
	Words are implemented as subclasses of :func:`tuple <py3:tuple>` meaning they are immutable.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	"""
	#Creation
	def __new__(cls, letters, arity, alphabet_size):
		"""Creates a new Word from a list of integers representing letters. The arity and alphabet size are stored. 
		The word's *letters* may be given as a list of integers or as a string (which is passed to the :func:from_string: function). Words are :func:`validated <validate>` and :func:`reduced to standard form <standardise>` before they are stored as a tuple.
		
		:raises TypeError: if *letters* is neither a string nor a list.
		:raises ValueError: see the errors raised by :func:`validate`.
		
			>>> Word("x2 a1 a2 a3 x1 a1 a2 a1 x2 L a2", 3, 2)
			Word('x1 a1 a2 a1', 3, 2)
		"""
		return cls._internal_new(letters, arity, alphabet_size)
	
	@classmethod
	def _internal_new(cls, letters, arity, alphabet_size, preprocess=True):
		"""Allows internal methods to bypass the validation checks."""
		if isinstance(letters, str):
			letters = from_string(letters)
		
		if preprocess:
			#1. Check to see that this has the right pattern of subwords and lambdas.
			#   Check that the letters are within the range -arity .... alphabet_size inclusive.
			validate(letters, arity, alphabet_size)
			#2. Standardise the list of symbols: remove as many lambdas as possible.
			letters = standardise(letters, arity)
		
		self = tuple.__new__(cls, letters)
		self.arity = arity
		self.alphabet_size = alphabet_size
		return self
	
	#Representation
	def __str__(self):
		return format(self)
	
	def __repr__(self):
		return "Word('{}', {}, {})".format(str(self), self.arity, self.alphabet_size)
	
	#Comparison
	def __lt__(self, other):
		r"""Words in standard form are compared according to dictionary order. The alphabet :math:`X \cup \Omega` is ordered according to:
		
		..math :: \dots < \alpha_2 < \alpha_1 < \lambda < x_1 < x_2 < \dots
		
			>>> Word('x1')
		"""
		#TODO examples and doctests
		if not isinstance(other, Word):
			return NotImplemented
		
		for s, o in zip(self, other):
			if s == o:
				continue
			if s < 0 and o < 0: #two alphas
				#self < o iff the index for s is less than the index for o
				return -s < -o
			return s < o
		return False
