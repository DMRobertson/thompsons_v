r"""By a *word*, we mean a string wrtitten using the symbols

.. math:: X \cup \Omega = \{x_1, \dotsc, x_r\} \cup \{\alpha_1, \dotsc, \alpha_n, \lambda\}.

The :math:`x_i` are constants, the :math:`\alpha_i` are unary operators and :math:`\lambda` is an *n*-ary operation. In this package, we refer to the parameters :math:`n` and :math:`r` as the *arity* and *alphabet size* respectively.

We consider three different algebras (in the sense of `universal algebra <http://en.wikipedia.org/wiki/Universal_algebra>`_) of words. In the notation of [Cohn]_:

1. :math:`W(\Omega; X)`, the set of all finite strings written using letters in :math:`X \cup \Omega`. Cohn calls these :math:`\Omega`-rows.
2. :math:`W_\Omega(X)`, the subset of :math:`W(\Omega; X)` whose strings begin with an :math:`x_i`, and represent a :func:`valid <validate>` series of operations. Cohn calls these :math:`\Omega`-words.
3. :math:`V_{n, r}(X) = W_\Omega(X) / \mathfrak{q}`. This is equivalent to the set of words in :math:`W_\Omega(X)` which are in Higman's :func:`standard form <standardise>`. Roughly speaking, a is in standard form if is valid (type 2) and it cannot be reduced to a shorter word using the following two rules. 

	a. :math:`w_1 \dots w_n \lambda \alpha_i = w_i`
	b. :math:`w \alpha_1 \dots w \alpha_n \lambda = w`

In the rules above, :math:`w` and the :math:`w_i` are all words in standard form.

Sometimes we need to refer to a string which consists only of :math:`\alpha`-s. Write :math:`A` for the set :math`\{\alpha_1, \dotsc, \alpha_n\}`. We define :math:`\langle A \rangle` to be the set of finite strings over :math:`A`.

.. seealso:: Definition 2.1, Section 2.3, Remark 3.3 and Definition 3.7 of the paper.

In the workings of this module, we represent a word as a :class:`tuple <py3:tuple>` of integers, where:

- :math:`x_i` is represented by :math:`i`,
- :math:`\alpha_i` is represented by :math:`-i`, and
- :math:`\lambda` is represented by :math:`0`.

We can write words of all types in this format, but we're only interested in the standard forms (type 3). To this end, the :func:`validate` detects those which are of type 2, and :func:`standardise` reduces type 2 words to type 3.

.. testsetup:: 
	
	from thompson.word import *

Functions and the Word class
----------------------------
"""

import operator
from itertools import chain
from .full_tree import FullTree

__all__ = ["format", "from_string", "validate", "standardise", "are_contractible", "lambda_arguments", "Word", "_concat"]

def _char(symbol):
	"""Converts the internal representation of a symbol into a nice string."""
	if symbol > 0:
		return 'x' + str(symbol)
	elif symbol < 0:
		return 'a' + str(-symbol)
	elif symbol == 0:
		return 'L'

def format(letters):
	"""Turns a sequence of integers representing a word into a nicely formatted string. Can be thought of as an inverse to :func:`from_string`. The word is not processed or reduced in any way.
	
		>>> format([2, -1, 2, -2, 0])
		'x2 a1 x2 a2 L'
	"""
	if len(letters) == 0:
		return "<the empty word>"
	return " ".join(_char(i) for i in letters)

def from_string(str):
	"""Converts a string representing word to the internal format---a tuple of integers. Anything which does not denote a basis letter (e.g. ``'x'``, ``'x2'``, ``'x45'``) a descendant operator (e.g.  ``'a1'``, ``'a3'``, ``'a27'``) or a  contraction (``'L'`` for :math:`\lambda`) is ignored.
	
	Note that ``'x'`` is interpreted as shorthand for ``'x1'``.
	
		>>> x = from_string('x2 a1 x2 a2 L')
		>>> print(x); print(format(x))
		(2, -1, 2, -2, 0)
		x2 a1 x2 a2 L
	
	:rtype: :class:`tuple <py3:tuple>` of integers
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
	return tuple(output)

def validate(letters, arity, alphabet_size=float('inf')):
	r"""Checks that the given *letters* describe a valid word. If no errors are raised when running this function, then we know that *letters* represents a valid word. Specifically, we know that:
	
	- The first letter is an :math:`x_i`.
	- Every :math:`\lambda` has *arity* arguments to its left.
	- If thge word contains a :math:`\lambda`, every subword is an argument of some :math:`\lambda`.
	- No :math:`\alpha_i` appears with :math:`i >` *arity*.
	- No :math:`x_i` occurs with  :math:`i >` *alphabet_size*.
	
	  - This last check occurs only when *alphabet_size* is given.
	
	Calling this function does not modfiy *letters*, and thus *letters* many not be in standard form (type 3).
	
		>>> validate(from_string('a1 x a2 L'), arity=2)
		Traceback (most recent call last):
			...
		ValueError: The first letter (a1) should be an x.
		>>> validate(from_string('x1 a2 x12 a1 a2 L'), arity=2, alphabet_size=4)
		Traceback (most recent call last):
			...
		IndexError: Letter x12 at index 2 is not in the alphabet (maximum is x4).
		>>> validate(from_string('x1 a1 a2 a3 a4 a5 x2 a2 L'), arity=3)
		Traceback (most recent call last):
			...
		IndexError: Letter a4 at index 4 is not in the alphabet (maximum is a3).
		>>> validate(from_string('x1 a2 L'), arity=4)
		Traceback (most recent call last):
			...
		ValueError: Letter L at index 2 is invalid. Check that lambda has 4 arguments.
		>>> validate(from_string('x1 x2 L x3 L x4'), arity=2)
		Traceback (most recent call last):
			...
		ValueError: Word is invalid: valency is 2 (should be 1).
	
	:raises ValueError: if this word fails the valency test.
	:raises IndexError: if this word contains an :math:`x_i` outside of the range 1 ... *alphabet_size*
	:raises IndexError: if this word contains an :math:`\alpha_i` outside of the range 1 ... *arity*
	:raises IndexError: if this word is empty (i.e. consists of 0 letters).
	
	.. seealso:: Proposition 2.12 of the paper for the *valency test*.
	"""
	symbol = letters[0]
	if symbol < 0:
		raise ValueError("The first letter ({}) should be an x.".format(
		  _char(symbol)))
	valency = 0
	for i, symbol in enumerate(letters):
		if symbol > alphabet_size:
			raise IndexError("Letter {} at index {} is not in the alphabet (maximum is x{}).".format(
			  _char(symbol), i, alphabet_size))
		if symbol < -arity:
			raise IndexError("Letter {} at index {} is not in the alphabet (maximum is a{}).".format(
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

def standardise(letters, arity, tail=()):
	"""Accepts a valid word as a :class:`tuple <py3:tuple>` of *letters* and reduces it to standard form. The result---a new tuple of letters---is returned. The *arity* must be given so we know how many arguments to pass to each :math:`\lambda`. The *tail* argument is used internally in recursive calls to this function and should be omitted.
	
	In the examples below, the :class:`Word` class standardises its input before creating a Word.
	
		>>> #Extraction only
		>>> print(Word("x a1 a2 a1 x a2 a1 L a1 a2", 2, 1))
		x1 a1 a2 a1 a2
		>>> #Contraction only
		>>> print(Word("x2 a2 a2 a1 x2 a2 a2 a2 L", 2, 2))
		x2 a2 a2
		>>> #Contraction only, arity 3
		>>> print(Word("x1 a1 a1 x1 a1 a2 x1 a1 a3 L", 3, 2))
		x1 a1
		>>> #Both extraction and contraction
		>>> print(Word("x a1 a2 a1 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L", 2, 1))
		x1 a1 a2 a1
		>>> #Lambda concatenation
		>>> print(Word("x a1 a2 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L", 2, 1))
		x1 a1 a2 a1 x1 a1 a2 a1 a2 L
		>>> #A mix of contraction and extraction
		>>> print(Word("x a2 x a1 L x a1 L a1 a2 a1 a2", 2, 1))
		x1 a1 a1 a2
		>>> #Can't contract different x-es
		>>> print(Word('x1 a1 x2 a2 L', 2, 2))
		x1 a1 x2 a2 L
		>>> #Something meaty, arity 3
		>>> print(Word('x1 a1 x1 a2 x1 a3 L a2 a1 a3 x1 a2 a1 x2 a1 x1 a1 a1 x1 a1 a2 x2 L x1 a2 L a2 a1 a1 L a3 a3 a2', 3, 2))
		x1 a1 a1 a1 a3 a2
	
	:raises TypeError: if *letters* is not a tuple of integers.
	:raises IndexError: if *letters* describes an :func:`invalid <validate>` word.
	
	:rtype: :class:`tuple <py3:tuple>` of integers.
	
	.. seealso:: Remark 3.3 of the paper.
	"""
	if not isinstance(letters, tuple):
		raise TypeError('Letters should be a tuple or Word instance. Instead, received {:r}'.format(
		  type(letters)))
	#0. Assume that *letters* describes a valid word, possibly not in standard form.
	#1. Find the rightmost lambda in the string.
	for i, symbol in enumerate(reversed(letters)):
		if symbol == 0:
			if i != 0:
				tail = letters[-i:] + tail
				letters = letters[:-i]
			break
	else:
		# If there is no such lambda, we have a simple word.
		return tuple(letters + tail)
	
	#2. Otherwise, break into the arguments of the lambda.
	subwords = lambda_arguments(letters, arity)
	
	#3. If the lambda is followed by alphas, extract subword and append any remaining alphas.
	if tail:
		alpha, tail = -tail[0], tail[1:]
		letters = subwords[alpha-1]
		return standardise(letters, arity, tail)
	
	#4. Otherwise, standardise each subword and attempt a lambda contraction.
	subwords = [standardise(subword, arity) for subword in subwords]
	prefix = are_contractible(subwords)
	if prefix:
		return prefix
	return _concat(subwords)

def lambda_arguments(word, arity=None):
	"""This function takes a :func:`valid <validate>` word which ends in a :math:`\lambda`, and returns the arguments of the rightmost :math:`\lambda` as a list of :class:`Words <Word>`.
	
	If *word* is given as a :class:`Word` instance, then *arity* may be omitted. Otherwise the *arity* should be provided.
	
		>>> w = 'x a1 a2 x a2 a2 L x a1 a1 L'
		>>> subwords = lambda_arguments(Word(w, 2, 1))
		>>> for subword in subwords: print(subword)
		x1 a1 a2 x1 a2 a2 L
		x1 a1 a1
		>>> w = 'x x x x a1 x L x a1 x a2 x L L x a1 a2 x L x a1 a2 x a2 a1 x L L'
		>>> subwords = lambda_arguments(Word(w, 3, 1))
		>>> for subword in subwords: print(format(subword))
		x1
		x1 x1 x1 a1 x1 L x1 a1 x1 a2 x1 L L x1 a1 a2 x1 L
		x1 a1 a2 x1 a2 a1 x1 L
	
	:raises IndexError: if *word* is an empty list of letters.
	:raises ValueError: if the last letter in *word* is not a lambda.
	:raises TypeError: if no arity is provided and *word* has no arity attribute.
	
	:rtype: list of :class:`Words <Word>`.
	"""
	is_Word_instance = isinstance(word, Word)
	if is_Word_instance and arity is None:
		arity = word.arity
	if word[-1] != 0:
		raise ValueError('The last letter of `...{}` is not a lambda.'.format(
		  format(word[-5:])))
	#Go backwards through the string. Track the largest partial sum of valencies.
	#When a new maxmium is found, you have just found the end of a subword.
	valency = 1 - arity
	max_valency = valency
	subword_end = -1
	subwords = []
	for i, symbol in enumerate(reversed(word[:-1])):
		valency += _valency_of(symbol, arity)
		if valency > max_valency:
			max_valency = valency
			subwords.append(word[-i - 2 : subword_end])
			subword_end = -i - 2
	
	assert i + 2 == len(word)
	assert len(subwords) == arity
	assert max_valency == valency == 1, (max_valency, valency)
	subwords.reverse()
	if is_Word_instance:
		for i in range(len(subwords)):
			subwords[i] = Word(subwords[i], word.arity, word.alphabet_size, preprocess=False)
	
	return subwords

def are_contractible(words):
	r"""Let *words* be a list of words, either as sequences of integers or as full :class:`Word` objects. This function tests to see if *words* is a list of the form :math:`(w\alpha_1, \dotsc, w\alpha_n)`, where ``n == len(words)``.
	
	:return: :math:`w` (as a :class:`tuple <py3:tuple>` of integers) if the test passes; the empty tuple if the test fails.
	
	**NB:** If there is a prefix :math:`w`, this function does **not** check to see if
	
	- :math:`w` :func:`is a valid word <validate>`, or
	- :math:`n` is the same as the arity of the context we're working in.
	
	>>> prefix = are_contractible([Word("x a1 a2 a3 a1", 3, 1), Word("x a1 a2 a3 a2", 3, 1), Word("x a1 a2 a3 a3", 3, 1)])
	>>> format(prefix)
	'x1 a1 a2 a3'
	"""
	#REMARK. Slicing a Word gives a tuple. This came about as a consequence of inheriting from tuple, but it makes sense.  Word objects are supposed to be fixed in standard form, and having a tuple instead of a word says that "this won't necessarily be a word".
	prefix = words[0][:-1]
	if len(prefix) == 0:
		return ()
	expected_length = len(prefix) + 1
	
	for j, word in enumerate(words):
		if not (
		  len(word) == expected_length
		  and word[:len(prefix)] == prefix
		  and word[-1] == -j - 1): #alpha_{j+1}
			return ()
	return prefix

def _concat(words):
	"""Takes an iterable *words* which yields sequences of integers representing words. Returns a tuple containing all the *words* concatenated together, with a zero (lambda) added on the end."""
	return tuple(chain.from_iterable(words)) + (0,)

class Word(tuple):
	r"""A sequence of letters (stored as integers), together with a recorded *arity* and *alphabet_size*. Words are implemented as subclasses of :class:`tuple <py3:tuple>`, meaning they are `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_. This means they can be used as dictionary keys. 
	
		>>> w = Word('x1 a1', 2, 1)
		>>> x = {}; x[w] = 'stored value'
		>>> x[w]
		'stored value'
		>>> #Can also use the underlying tuple as a key
		>>> x[1, -1]
		'stored value'
		>>> #the reason this works
		>>> hash(w) == hash((1, -1))
		True
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	:ivar lambda_length: the number of times the contraction symbol :math:`\lambda` occurs in this word after it has been written in standard form.
	"""
	
	#I'd previously thought about defining __slots__ to save memory. But Python doesn't like that on subclasses of builtin or immutable types (I forget which).
	#Creation
	def __new__(cls, letters, arity, alphabet_size, preprocess=True):
		"""Creates a new Word consisting of the given *letters* and stores the *arity* and *alphabet_size*. 
		The *letters* may be given as a list of integers or as a string (which is passed to the :func:`from_string` function).
		
			>>> Word("x2 a1 a2 a3 x1 a1 a2 a1 x2 L a2", 3, 2)
			Word('x1 a1 a2 a1', 3, 2)
		
		By default, ``preprocess is True``. This means that words are :func:`validated <validate>` and :func:`reduced to standard form <standardise>` before they are stored as a tuple. These steps can be disabled by using ``preprocess=False``. This option is used internally when we *know* we have letters which are already vaild and in standard form.
		
		:raises ValueError: See the errors raised by :func:`validate`.
		"""
		if isinstance(letters, str):
			letters = from_string(letters)
		elif not isinstance(letters, tuple):
			letters = tuple(letters)
		
		if preprocess:
			#1. Check to see that this has the right pattern of subwords and lambdas.
			#   Check that the letters are within the range -arity .... alphabet_size inclusive.
			validate(letters, arity, alphabet_size)
			#2. Standardise the letters: remove as many lambdas as possible.
			letters = standardise(letters, arity)
		
		self = tuple.__new__(cls, letters)
		self.arity = arity
		self.alphabet_size = alphabet_size
		self.lambda_length = self.count(0)
		return self
	
	#Representation
	def __str__(self):
		return format(self)
	
	def __repr__(self):
		return "Word('{}', {}, {})".format(str(self), self.arity, self.alphabet_size)
	
	#Comparisons
	def __lt__(self, other):
		r"""Let us assign a total order to the alphabet :math:`X \cup \Omega` by declaring:
		
		.. math:: \alpha_1 < \alpha_2 < \alpha_3 < \dots < x_1 < x_2 < x_3 < \dots < \lambda 
		
		We can use this to order :meth:`simple words <is_simple>` by using `dictionary order <http://en.wikipedia.org/wiki/Lexicographical_order#Motivation_and_uses>`_.
		
		>>> Word('x1', 2, 2) < Word('x1', 2, 2)
		False
		>>> Word('x1', 2, 2) < Word('x2', 2, 2)
		True
		>>> Word('x1 a2', 2, 2) < Word('x1 a1', 2, 2)
		False
		>>> Word('x1 a1 a2 a3 a4', 4, 1) < Word('x1 a1 a2 a3 a3', 4, 1)
		False
		
		We extend this to words :func:`in standard form <standardise>` in the following way. Let :math:`\lambda(u)` denote the lambda-length of :math:`u`.
		
		1. If :math:`\lambda(u) < \lambda(v)`, then :math:`u < v`.
		
		>>> Word('x2 x1 L', 2, 2) < Word('x1 x2 x1 L L', 2, 2)
		True
		
		2. If :math:`u=v`, then neither is strictly less than the other.
		
		>>> Word('x1 x2 L', 2, 2) < Word('x1 x2 L', 2, 2)
		False
		
		3. Otherwise :math:`u \neq v` are different words with the same :math:`\lambda`-length. Break both words into the :math:`n` arguments of the outmost lambda, so that :math:`u = u_1 \dots u_n \lambda` and :math:`v = v_1 \dots v_n \lambda`, where each subword is in standard form.
		4. Let :math:`i` be the first index for which :math:`u_i \neq v_i`.  Test if :math:`u_i < v_i` by applying this definition recursively. If this is the case, then :math:`u < v`; otherwise, :math:`u > v`.
		
		>>> #True, because x2 < x2 a2
		>>> Word('x1 x2 L', 2, 2) < Word('x1 x2 a2 L', 2, 2)
		True
		>>> #True, because words of lambda-length 1 are less than those of lambda-length 2.
		>>> Word('x1 x2 L x3 x4 L L', 2, 4) < Word('x1 x2 L x3 L x4 L', 2, 4)
		True
		
		The other three comparison operators (``<=, >, >=``) are also implemented.
		
		>>> Word('x1 x2 L', 2, 2) <= Word('x1 x2 L', 2, 2)
		True
		>>> Word('x1 a1', 2, 2) <= Word('x1 a1', 2, 2) <= Word('x1 a1 a2', 2, 2)
		True
		>>> Word('x1', 2, 2) > Word('x1', 2, 2)
		False
		>>> Word('x1', 2, 2) >= Word('x1', 2, 2)
		True
		>>> Word('x1 a2', 2, 2) > Word('x1', 2, 2)
		True
		>>> Word('x1 a2', 2, 2) >= Word('x1', 2, 2)
		True
		
		.. todo:: I think this is a total order---try to prove this. I based the idea on [Zaks]_ (section 2, definition 1) which describes a total order on *k*-ary trees.
		"""
		return self._compare(other, operator.lt, False)
	
	def __gt__(self, other):
		return self._compare(other, operator.gt, False)
	
	def __le__(self, other):
		return self._compare(other, operator.lt, True)
	
	def __ge__(self, other):
		return self._compare(other, operator.gt, True)
	
	def _compare(self, other, comparison, nonstrict):
		"""The recursive part of the less than/greater than test. See the docstring for __lt__.
		
		- Comparison should be one of operator.lt and operator.gt.
		- For strict inequality, pass nonstrict = False
		- For nonstrict inequality, pass nonstrict = True
		"""
		if not isinstance(other, Word):
			return NotImplemented
		#Simple words are less than words with 1 lambda; they are less than words with 2 lambdas, etc.
		if self.lambda_length != other.lambda_length:
			return comparison(self.lambda_length, other.lambda_length)
		
		#Simple words are done with dictionary order.
		if self.lambda_length == 0:
			return self._compare_simple(other, comparison, nonstrict)
		
		#otherwise we have two words 
		self_args  = lambda_arguments(self)
		other_args = lambda_arguments(other)
		for s, o in zip(self_args, other_args):
			if s == o:
				continue
			else:
				return s._compare(o, comparison, nonstrict)
		return nonstrict
	
	def _compare_simple(self, other, comparison, nonstrict):
		"""The basis for the induction of the less than/greater than test. See the docstring for __lt__.
		
		>>> Word('x1', 2, 2) < Word('x2', 2, 2)
		True
		>>> Word('x1 a1', 2, 2) < Word('x1', 2, 2)
		False
		>>> Word('x1 a1', 2, 2) < Word('x1 a2', 2, 2)
		True
		"""
		assert self.lambda_length == other.lambda_length == 0, "_compare_simple called on non-simple arguments"
		for s, o in zip(self, other):
			if s == o:
				continue
			#Else s and o are different and not lambdas
			if s < 0 and o < 0: #two alphas
				return comparison(-s, -o)
			return comparison(s, o)
		
		#If we reach here, then both words start with the same thing.
		if len(self) == len(other):
			return nonstrict
		#Otherwise, one is longer than the other.
		return comparison(len(self), len(other))
	
	#Tests
	def is_simple(self):
		"""Let us call a Word *simple* if it has a lambda length of zero once it has been reduced to standard form.
		
			>>> Word("x1 a1 a2", 2, 1).is_simple()
			True
			>>> Word("x1 a1 a2 x2 a2 L", 2, 2).is_simple()
			False
			>>> #Simplify a contraction
			>>> Word("x1 a1 a1 x1 a1 a2 L", 2, 1).is_simple()
			True
			>>> #Lambda-alpha extraction
			>>> Word("x1 x2 L a2", 2, 2).is_simple()
			True
		"""
		return self.lambda_length == 0
	
	def is_above(self, word):
		""":meth:`Tests <test_above>` to see if the current word is an initial segment of *word*. Returns True if the test passes and False if the test fails.
		
			>>> w = Word('x1 a1', 2, 2)
			>>> w.is_above(Word('x1 a1 a1 a2', 2, 2))
			True
			>>> w.is_above(Word('x1', 2, 2))
			False
			>>> w.is_above(w)
			True
			>>> v = Word('x1 a1 x1 a2 L', 2, 2)
			>>> w.is_above(v)
			False
			>>> v.is_above(w)
			True
		"""
		return self.test_above(word) is not None
	
	def test_above(self, word):
		r"""Tests to see if the current word :math:`c` is an initial segment of the given *word* :math:`w`. In symbols, we're testing if :math:`w = c \Gamma`, where :math:`\Gamma \in \langle A \rangle` is some string of alphas only.
		
		The test returns :math:`\Gamma` (as a tuple of integers) if the test passes; note that :math:`\Gamma` could be empty (if :math:`c = w`). If the test fails, returns ``None``.
		
			>>> c = Word('x1 a1', 2, 2)
			>>> c.test_above(Word('x1 a1 a1 a2', 2, 2))
			(-1, -2)
			>>> c.test_above(Word('x1', 2, 2)) is None
			True
			>>> c.test_above(c)
			()
			>>> w = Word('x1 a2 x1 a1 L', 2, 2)
			>>> print(c.test_above(w))
			None
			>>> w.test_above(c)
			(-2,)
			>>> v = Word('x1 a2 x1 a1 L x1 a2 a2 L x1 a2 x1 a1 L L', 2, 2)
			>>> print(c.test_above(v))
			None
			>>> #There are two possible \Gamma values here; only one is returned.
			>>> v.test_above(c)
			(-1, -1, -2)
			
		.. note:: There may be many different values of :math:`\Gamma` for which :math:`w = c \Gamma`; this method simply returns the first such :math:`\Gamma` that it finds.
		"""
		if self.lambda_length < word.lambda_length:
			return None
		
		if self.lambda_length == word.lambda_length > 0:
			#self above word iff self == word
			return tuple() if self == word else None
		
		if self.lambda_length == word.lambda_length == 0:
			if len(self) >  len(word):
				return None
			if len(self) <= len(word):
				head, tail = word.split(len(self))
				return tail if head == self else None
		
		if self.lambda_length > word.lambda_length:
			subwords = lambda_arguments(self)
			for i, subword in enumerate(subwords):
				tail = subword.test_above(word)
				if tail is not None:
					return (-(i + 1),) + tail
			return None
	
	def max_depth_to(self, basis):
		r"""Let :math:`w` be the current word and let :math:`X` be a :meth:`basis <thompson.generators.Generators.is_basis>`. Choose a (possibly empty) string of alphas :math:`\Gamma` of length :math:`s \geq 0` at random. What is the smallest value of :math:`s` for which we can guarantee that :math:`w\Gamma` is below :math:`X`, i.e. in :math:`X\langle A\rangle`?
		
			>>> from thompson.generators import Generators
			>>> basis = Generators.standard_basis(2, 1).expand(0).expand(0).expand(0)
			>>> basis
			Generators(2, 1, ['x1 a1 a1 a1', 'x1 a1 a1 a2', 'x1 a1 a2', 'x1 a2'])
			>>> Word('x', 2, 1).max_depth_to(basis)
			3
			>>> Word('x a1', 2, 1).max_depth_to(basis)
			2
			>>> Word('x a2 x a1 L x1 a1 L', 2, 1).max_depth_to(basis)
			4
			>>> Word('x a2', 2, 1).max_depth_to(basis)
			0
		"""
		max_depth = 0
		dict = {self: 0}
		while dict:
			word, expansions = dict.popitem()
			if basis.is_above(word):
				max_depth = max(max_depth, expansions)
				continue
			else:
				for child in word.expand():
					dict[child] = expansions + 1
		return max_depth
	
	#Modifiers
	def alpha(self, index):
		r"""Let :math:`w` stand for the current word. This method creates and returns a new word :math:`w\alpha_\text{index}`.
		
			>>> Word("x a1 a2", 3, 2).alpha(3)
			Word('x1 a1 a2 a3', 3, 2)
		
		:raises IndexError: if *index* is not in the range 1... *arity*.
		"""
		if not 1 <= index <= self.arity:
			raise IndexError("Index ({}) is not in the range 1...{}".format(
			  index, self.arity))
		preprocess = not self.is_simple()
		return type(self)(self + (-index,), self.arity, self.alphabet_size, preprocess)
	
	def extend(self, tail):
		"""Concatenates the current word with the series of letters *tail* to form a new word.The argument *tail* can be given as either a string or a tuple of integers.
		
			>>> Word('x1 a2 a1', 2, 1).extend('a2 a2')
			Word('x1 a2 a1 a2 a2', 2, 1)
			>>> Word('x1 a2 a1', 2, 1).extend('x1 a2 a2')
			Traceback (most recent call last):
				...
			ValueError: Word is invalid: valency is 2 (should be 1).
			>>> Word('x1 a2 a1', 2, 1).extend('x1 a2 a2 L')
			Word('x1 a2', 2, 1)
		"""
		if isinstance(tail, str):
			tail = from_string(tail)
		return type(self)(self + tail, self.arity, self.alphabet_size)
	
	def expand(self):
		"""Returns an iterator that yields the *arity* descendants of this word.
		
			>>> w = Word("x a1", 5, 1)
			>>> for child in w.expand():
			... 	print(child)
			x1 a1 a1
			x1 a1 a2
			x1 a1 a3
			x1 a1 a4
			x1 a1 a5
			>>> w = Word("x a1 a1 x a2 a1 L", 2, 1)
			>>> for child in w.expand():
			... 	print(child)
			x1 a1 a1
			x1 a2 a1
		
		:rtype: iterator which yields :class:`Words <Word>`.
		"""
		if self.is_simple():
			return (self.alpha(i) for i in range(1, self.arity + 1))
		return (type(self)(subword, self.arity, self.alphabet_size, preprocess=False)
		  for subword in lambda_arguments(self))
	
	def _as_contraction(self):
		r"""Converts a Word :math:`w` in standard form to the equivalent word :math:`w \alpha_1 \dots w \alpha_n \lambda`. The new word **is not** in standard form, though it is valid. This is *reaaaaaalllly naughty*, because Word instances are meant to be always in standard form. But it has some uses, honest. Just be careful with this one.
		"""
		letters = []
		for i in range(self.arity):
			letters.extend(self)
			letters.append(-(i + 1))
		letters.append(0)
		return type(self)(letters, self.arity, self.alphabet_size, preprocess=False)
	
	def split(self, head_width):
		"""Splits the current word *w* into a pair of tuples *head*, *tail* where ``len(head) == head_width``. The segments of the word are returned as tuples of integers, i.e. not fully fledged :class:`Words <Word>`.
		
		:raises IndexError: if *head_width* is outside of the range ``0 <= head_width <= len(w)``.
		
		>>> Word('x1 a1 a2 a3 a1 a2', 3, 1).split(4)
		((1, -1, -2, -3), (-1, -2))
		>>> Word('x1 a1 a2 a3 a1 a2', 3, 1).split(10)
		Traceback (most recent call last):
			...
		IndexError: The head width 10 is not in the range 0 to 6.
		"""
		if not 0 <= head_width <= len(self):
			raise IndexError('The head width {} is not in the range 0 to {}'.format(
			  head_width, len(self)))
		return self[:head_width], self[head_width:]
	
	def rsplit(self, tail_width):
		"""The same as :meth:`split`, but this time *tail_width* counts from the right hand side. This means that ``len(tail) == tail_width``.
		
		:raises IndexError: if *tail_width* is outside of the range ``0 <= tail_width <= len(w)``.
		
		>>> Word('x1 a1 a2 a3 a1 a2', 3, 1).rsplit(4)
		((1, -1), (-2, -3, -1, -2))
		>>> Word('x1 a1 a2 a3 a1 a2', 3, 1).rsplit(10)
		Traceback (most recent call last):
			...
		IndexError: The tail width 10 is not in the range 0 to 6.
		"""
		if not 0 <= tail_width <= len(self):
			raise IndexError('The tail width {} is not in the range 0 to {}'.format(
			  tail_width, len(self)))
		return self[:-tail_width], self[-tail_width:]

