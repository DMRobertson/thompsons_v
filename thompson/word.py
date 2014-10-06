r"""The algebra :math:`V_{n, r}(X)` consists of words written in the alphabet :math:`X \cup \Omega = \{x_1, \dotsc, x_r\} \cup \{\alpha_1, \dotsc, \alpha_n, \lambda\}`, with certain restrictions (defined in remark 3.3). In fact, there are three different algebras of words. Using the notation of [Cohn]_:

	1. :math:`W(\Omega; X)`, the set of finite strings written using letters in  :math:`X \cup \Omega`. Cohn calls these :math:`\Omega` *-rows*.
	2. :math:`W_\Omega(X)`, the subset of :math:`W(\Omega; X)` whose strings begin with an :math:`x_i`, and represent a :func:`valid <validate>` series of operations. Cohn calls these :math:`\Omega` *-words*.
	3. :math:`V_{n, r}(X) = W_\Omega(X) / \mathfrak{q}`. This is (equivalent to) the set of words in :math:`W_\Omega(X)` which are in Higman's :func:`standard form <standardise>` (see remark 3.3 and lemma 3.4).

In the workings of this module, we represent a word as a :class:`tuple <py3:tuple>` (or any type of sequence) of integers, where:

- :math:`x_i` is represented by :math:`i`,
- :math:`\alpha_i` is represented by :math:`-i`, and
- :math:`\lambda` is represented by :math:`0`.

We can write words of all types in this format, but we're only interested in the standard forms (type 3). To this end, the :func:`validate` detects those which are of type 2, and :func:`standardise` reduces type 2 words to type 3.

.. testsetup:: 
	
	from thompson.word import *
"""

import operator
from itertools import chain
from .full_tree import FullTree

__all__ = ["format", "from_string", "validate", "standardise", "are_contractible", "lambda_arguments", "Word", "initial_segment"]

def _char(symbol):
	if symbol > 0:
		return 'x' + str(symbol)
	elif symbol < 0:
		return 'a' + str(-symbol)
	elif symbol == 0:
		return 'L'

def format(letters):
	"""Turns a list of integers representing a word into a nicely formatted string. Can be thought of as an inverse to :func:`from_string`. No processing is applied to the word.
	
		>>> format([2, -1, 2, -2, 0])
		'x2 a1 x2 a2 L'
	"""
	return " ".join(_char(i) for i in letters)

def from_string(str):
	"""Converts a string representing word to the internal format---a tuple of integers. Anything which does not denote a basis letter (e.g. ``'x'``, ``'x2'``, ``'x45'``) a descendant operator (e.g.  ``'a1'``, ``'x3'``, ``'x27'``) or a lambda contraction (``'L'``) is ignored.
	
		>>> x = from_string('x2 a1 x2 a2 L')
		>>> print(x); print(format(x))
		(2, -1, 2, -2, 0)
		x2 a1 x2 a2 L
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
	r"""Checks that the given *letters* describes a valid word (type 2), using the *valency test* of proposition 2.12 (due to [Cohn]_). If no errors are raised when running this function, then we know that *letters* represents a valid word. Specifically, we know that:
	
	- The first letter is an :math:`x_i`.
	- Every :math:`\lambda`-contraction has *arity* arguments to its left.
	- No :math:`\alpha_i` appears with :math:`i >` *arity*.
	- No :math:`x_i` occurs with  :math:`i >` *alphabet_size*. *NB:* this last check occurs only when *alphabet_size* is given.
	
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
	
	:raises ValueError: if this word fails the valency test (prop 2.12).
	:raises IndexError: if this word contains an :math:`x_i` outside of the range 1 ... *alphabet_size*
	:raises IndexError: if this word contains an :math:`\alpha_i` outside of the range 1 ... *arity*
	:raises IndexError: if this word is empty (i.e. consists of 0 letters).
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

def standardise(letters, arity):
	"""Accepts a word (as a list of letters) and reduces it to standard form. The result (a new list of letters) is then returned. See also remark 3.3.
	
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
	
	:raises IndexError: if the list of letters :func:`is not valid <validate>`.
	:raises TypeError: if *letters* is not a list of integers or a Word.
	"""
	#Form the tree of lambda calls
	root = _syntax_tree(letters, arity)
	
	#Modify the tree. Find all lambda-alpha nodes and replace them with their extracted child.
	extractions = []
	for node in root.walk_postorder():
		if not node.is_leaf() and node.letters:
			extractions.append(node)
	
	while extractions:
		node = extractions.pop()
		head, tail = node.letters[0], node.letters[1:]
		
		#The extracted child
		replacement = node.children[-head - 1]
		#Remove the discarded children from the extractions list if they're in there.
		for child in node:
			if child is not replacement:
				_remove_with_descendants(child, extractions)
		
		node.replace_with_child(-head - 1)
		if replacement.is_root():
			root = replacement 
		replacement.letters += tail
		
		if (replacement not in extractions
		  and replacement.letters 
		  and not replacement.is_leaf()):
			extractions.append(replacement)
	
	#Now the only lambdas left should be contractions.
	assert all(not node.letters for node in root.walk() if not node.is_leaf()), 'branch has letters'
	
	for node in root.walk_postorder():
		if node.is_leaf():
			continue
		subwords = [child.letters for child in node]
		prefix = are_contractible(subwords)
		if prefix: #is not empty, so can contract
			node.letters = prefix
		else: #we'll just have to concatenate won't we?
			node.letters = concat(subwords) + [0]
	
	letters = root.letters
	root.discard()
	return tuple(letters)

def _syntax_tree(letters, arity):
	"""Returns a tree representing the composition of lambdas in this word."""
	alphas = []
	indices = [0] #the indices of the rightmost unconsidered node on each level.
	root = node = FullTree(arity)
	
	#Break eveything down into a tree of lambda calls.
	#The leaves are the simple words in x<A>.
	for symbol in reversed(letters):
		if symbol < 0: #alphas
			alphas.append(symbol)
		
		elif symbol == 0: #lambda
			node.letters = alphas[::-1]
			node.expand()
			indices[-1] -= 1
			alphas = []
			indices.append(arity - 1)
			node = node.children[indices[-1]]
		
		elif symbol > 0: #x
			node.letters = [symbol] + alphas[::-1]
			indices[-1] -= 1
			alphas = []
			
			#Go up the tree if needed
			while indices and indices[-1] < 0:
				node = node.parent
				indices.pop()
			if indices and indices[-1] >= 0:
				#move to the left sibling
				node = node.parent.children[indices[-1]]
			
	assert node is None
	return root

def _remove_with_descendants(node, list):
	for descendant in node.walk_postorder():
		try:
			list.remove(descendant)
		except ValueError:
			pass

def are_contractible(words):
	r"""Let *words* be a list of words, either as lists of integers or as full :class:`Word` objects. This function tests to see if *words* is a list of the form :math:`(w\alpha_1, \dotsc, w\alpha_n)`, where ``n == len(words)``.
	
	:return: :math:`w` (as a tuple of integers) if the test passes; the empty tuple if the test fails.
	
	**NB:** If there is a prefix :math:`w`, this function does **not** check to see if :math:`w` :func:`is a valid word <validate>`.
	
		>>> prefix = are_contractible([Word("x a1 a2 a3 a1", 3, 1), Word("x a1 a2 a3 a2", 3, 1), Word("x a1 a2 a3 a3", 3, 1)])
		>>> format(prefix)
		'x1 a1 a2 a3'
	"""
	#REMARK. Slicing a Word gives a tuple. This came about as a consequence of inheriting from tuple, but it makes sense.  Word objects are supposed to be fixed in standard form, and having a tuple instead of a word says that "this won't neccesarily be a word".
	prefix = words[0][:-1]
	if len(prefix) == 0:
		return []
	expected_length = len(prefix) + 1
	
	for j, word in enumerate(words):
		if not (
		  len(word) == expected_length
		  and word[:len(prefix)] == prefix
		  and word[-1] == -j - 1): #alpha_{j+1}
			return []
	return prefix

def concat(iterable):
	return list(chain.from_iterable(iterable))

class Word(tuple):
	r"""A tuple of letters (stored as integers), together with a recorded *arity* and *alphabet_size*. Words are implemented as subclasses of :class:`tuple <py3:tuple>`, meaning they are `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_. This means they can be used as dictionary keys. 
	
		>>> w = Word('x1 a1', 2, 1)
		>>> x = {}; x[w] = 'stored value'
		>>> x[w]
		'stored value'
		>>> #Try using the underlying tuple as a key
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
		"""Creates a new Word from a sequence of letters and stores the *arity* and *alphabet_size*. 
		The word's *letters* may be given as a list of integers or as a string (which is passed to the :func:`from_string` function).
		
			>>> Word("x2 a1 a2 a3 x1 a1 a2 a1 x2 L a2", 3, 2)
			Word('x1 a1 a2 a1', 3, 2)
		
		By default, the argument *preprocess* is True. THis means that words are :func:`validated <validate>` and :func:`reduced to standard form <standardise>` before they are stored as a tuple. These steps can be disabled by using ``preprocess=False``; this is used internally when we know we have letters which are already vaild and in standard form.
		
		:raises TypeError: if *letters* is neither a string nor a list.
		:raises ValueError: See the errors raised by :func:`validate`.
		"""
		if isinstance(letters, str):
			letters = from_string(letters)
		
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
		r"""Let us assign a total order to the alphabet :math:`X \cup \Omega` by:
		
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
		
		We extend this to words :func:`in standard form <standardise>` involving lambda-contractions in the following way. Let :math:`\lambda(u)` denote the lambda-length of :math:`u`.
		
		1. :math:`\lambda(u) < \lambda(v) \implies u < v`.
		
		>>> Word('x2 x1 L', 2, 2) < Word('x1 x2 x1 L L', 2, 2)
		True
		
		2. If :math:`u` and :math:`v` are equal, then neither is strictly less than the other.
		
		>>> Word('x1 x2 L', 2, 2) < Word('x1 x2 L', 2, 2)
		False
		
		3. Else, :math:`\lambda(u) \neq \lambda(v)`. Break :math:`u` and :math:`v` into the :math:`n` arguments of the outmost lambda, so that :math:`u = u_1 \dots u_n \lambda` and :math:`v = v_1 \dots v_n \lambda`, where each subword is in standard form. Let :math:`i` be the first index for which :math:`\lambda(u_i) \neq \lambda(v_i)`. 
		4. Test if :math:`\lambda(u_i) < \lambda(v_i)` by applying this definition recursively. If this is the case, then :math:`\lambda(u) < \lambda(v)`; otherwise, :math:`\lambda(u) > \lambda(v)`.
		
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
		
		.. todo:: I think this is a total order. I based the idea on [Zaks]_ (section 2, definition 1) which describes a total order on *k*-ary trees. I should try to prove that this is (hopefully) a total order on *labelled* *k*-ary trees.
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
		"""Tests to see if the current word is an initial segment of *word*. Returns True if the test passes and False if the test fails.
		
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
		r"""Tests to see if the current word :math:`c` is an initial segment of the given word :math:`w`. In symbols, we're testing if :math:`w = c \Gamma`, where :math:`\Gamma \in \langle A \rangle` is some string of alphas only.
		
		The test returns :math:`\Gamma` (as a tuple of integers) if the test passes; note that :math:`\Gamma` could be empty (if :math:`c = s`). If the test fails, returns ``None``.
		
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
		return type(self)(self + (-index,), self.arity, self.alphabet_size, preprocess=False)
	
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
		"""
		return (self.alpha(i) for i in range(1, self.arity + 1))
	
	def split(self, head_width):
		"""Splits the current word *w* into a pair of tuples *head*, *tail* where ``len(head) == head_width``. The segments of the word are returned as tuples of integers (not fully fledged words).
		
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

#Operations on words in standard form
def initial_segment(u, v):
	r"""Let :math:`u, v` be two words in standard form. This function returns True if :math:`u` is an initial segment of :math:`v` or vice versa; otherwise returns False. See definition 3.15.
	
	>>> #completely different words
	>>> u = Word("x1 a2 a1 a1 a2", 2, 2)
	>>> v = Word("x2 a1 a2 a1 a2", 2, 2)
	>>> initial_segment(u, v)
	False
	>>> #v starts with u
	>>> u = Word("x a1 a1 x2 a1 a2 L a1 a2", 2, 2)
	>>> v = Word("x a1 a1 a2 a2 a1", 2, 2)
	>>> initial_segment(u, v)
	True
	>>> #Before reduction, v starts with u. A later lambda changes this.
	>>> u = Word("x a1 a2", 2, 2)
	>>> v = Word("x a1 a2 x a2 a2 L", 2, 2)
	>>> initial_segment(u, v)
	False
	
	.. seealso:: :meth:`Word.is_above` does the same thing without the 'vice versa' part.
	"""
	if len(u) > len(v):
		u, v = v, u
	return all(a == b for a, b in zip(u, v)) and all(b < 0 for b in v[len(u):])

def lambda_arguments(word):
	"""This function takes a :class:`Word` (in standard form) which ends in a lambda, and returns the arguments of the rightmost lambda as a list of Words.
	
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
	
	:raises ValueError: if the last letter in *word* is not a lambda.
	:raises TypeError: if *word* is not a :class:`Word` instance.
	"""
	if not isinstance(word, Word):
		raise TypeError("The argument {} is not a Word instance".format(repr(word)))
	if word[-1] != 0:
		raise ValueError('The last letter of `...{}` is not a lambda.'.format(
		  format(word[-5:])))
	#Go backwards through the string. Track the largest partial sum of valencies.
	#When a new maxmium is found, you have just found the end of a subword.
	valency = 1 - word.arity
	max_valency = valency
	subword_end = -1
	subwords = []
	for i, symbol in enumerate(reversed(word[:-1])):
		valency += _valency_of(symbol, word.arity)
		if valency > max_valency:
			max_valency = valency
			subwords.append(word[-i - 2 : subword_end])
			subword_end = -i - 2
	
	assert i + 2 == len(word)
	assert len(subwords) == word.arity
	assert max_valency == valency == 1, (max_valency, valency)
	subwords.reverse()
	for i in range(len(subwords)):
		subwords[i] = Word(subwords[i], word.arity, word.alphabet_size, preprocess=False)
	
	return subwords




