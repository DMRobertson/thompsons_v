"""
.. testsetup:: 
	
	from thompson.word import *
"""

from itertools import chain

class Word(list):
	r"""Represents an element of the algebra :math:`V_{n, r}`. This consists of words written in the alphabet :math:`\{x_1,\dotsc, x_r, \alpha_1, \dotsc, \alpha_n, \lambda\}`. See Remark 3.3.
	"""
	#TODO is this really needed?
	def __init__(self, arity, alphabet_size):
		r"""When creating a word manually, you must specify the arity *n* and the number of :math:`x_i`, *r*."""
		super().__init__(self)
		
		self.arity = arity
		self.alphabet_size = alphabet_size
	
	#Internal scheme: positive integers  i correspond to x_i
	#                 negative integers -i correspond to alpha_i
	#                 zero corresponds to lambda
	
	def letter(self, index=1):
		r"""Writes the letter :math:`x_\text{index}` at the end of this word.
		
		:raises ValueError: if index is not in the range 1, ..., *alphabet_size*."""
		if not 1 <= index <= self.alphabet_size:
			raise ValueError("Index {} is not in the range 1 to alphabet size {}".format(index, self.alphabet_size))
		self.append(index)
		return self
	
	# x = letter
	
	def descend(self, index, power=1):
		r"""Writes the descend operator :math:`\alpha_\text{index}` at the end of this word *power* times.
		
		:raises ValueError: if index is not in the range 1, ..., *arity*."""
		if not 1 <= index <= self.arity:
			raise ValueError("Index {} is in the range 1 to arity {}".format(index, self.arity))
		for _ in range(power):
			self.append(-index)
		return self
	
	# alpha = descend
	
	def contract(self):
		r"""Writes the contraction operator :math:`\lambda` at the end of this word."""
		self.append(0)
		return self
	
	# lambda_ = contract
	
	def standardise(self):
		r"""Reduces this word to its standard form as defined in Rem 3.3.
		
		:raises IndexError: if the first character is not a letter :math:`x_i`.
		:raises ValueError: if a :math:`\lambda` is found with insufficient arguments to its left.
		
			>>> w = Word.from_string("x2 a1 a2 a1 x2 a1 a2 a2 L a2 a1")
			>>> w.standardise()
			>>> print(w)
			x2 a1 a2 a2 a1
			>>> w = Word.from_string("x a1 a2 a1 x a2 a1 L a1 a2")
			>>> w.standardise()
			>>> print(w)
			x a1 a2 a1 a2
			>>> w = Word.from_string("x a1 a2 a1 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L")
			>>> w.standardise()
			>>> print(w)
			x a1 a2 a1
			>>> w = Word.from_string("x a1 a2 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L")
			>>> w.standardise()
			>>> print(w)
			x a1 a2 a1 x a1 a2 a1 a2 L
		"""
		#TODO: an internal _is_standardised flag
		subwords = []
		i = 0
		while i < len(self):
			symbol = self[i]
			# print("{:2} {:2} {}".format(symbol, i, subwords))
			if symbol > 0: #Letter
				subwords.append([symbol])
			
			elif symbol == 0: #Lambda
				if len(subwords) < self.arity:
					raise ValueError("Lambda contraction expected {} args but received only {}.".format(
						self.arity, len(subwords)
					))
				#Peek ahead: is the lambda immediately followed by an alpha?
				next = self[i+1] if i + 1 < len(self) else 0
				if next < 0:
					# print('extract')
					i = self._extract(subwords, i, -next)
				else:
					# print('contract')
					i = self._contract(subwords, i)
			
			else: #Alpha
				subwords[-1].append(symbol)
				
			i += 1
	
	def _extract(self, subwords, lambda_position, next):
		"""If we find a lambda followed by an alpha, replace the lambda, the alpha and the lambda's arguments by argument #*next*.
		
		:return: Returns one less than what the value of i should be on the next iteration."""
		start = lambda_position - sum(len(list) for list in subwords[-self.arity:])
		extracted = subwords[-self.arity + next - 1]
		
		del subwords[-self.arity + 1:]
		subwords[-1] = extracted
		self[start:lambda_position + 2] = extracted
		return start + len(extracted) - 1
	
	@classmethod
	def are_contractible(cls, words):
		r"""Tests to see if *words* is a list of the form :math:`(w\alpha_1, \dotsc, w\alpha_n)`.
		
		:return: :math:`w` if the test passes; ``None`` if the test fails.
		"""
		prefix = words[0][:-1]
		assert prefix[0] > 0  #is a letter
		expected_length = len(prefix) + 1
		
		for j, word in enumerate(words):
			if not (
			  len(word) == expected_length
			  and word[:len(prefix)] == prefix
			  and word[-1] == -j - 1): #alpha_{j+1}
				return None
		return prefix
	
	def _contract(self, subwords, lambda_position):
		"""Attempts to contract the last *n* words in *subwords*, where *n* is the current word's arity.
		
		Returns one less than what the value of i should be on the next iteration.
		"""
		inspected = subwords[-self.arity:]
		start = lambda_position - sum(len(list) for list in inspected)
		prefix = Word.are_contractible(inspected)
		
		if prefix is not None:
			#TODO: The lambda call is not of the form ua1, ua2, ..., uan, lambda, and the next symbol is NOT an alpha.
			#I think there's no way this can be reduced further.
			subwords[-self.arity] = prefix
			del subwords[-self.arity + 1:]
			self[start : lambda_position + 1] = prefix
			return start + len(prefix) - 1
		else:
			#We can contract the last *arity* words.
			subwords[-self.arity] = chain.from_iterable(inspected)
			del subwords[-self.arity + 1:]
			return lambda_position
	
	@classmethod
	def from_string(cls, str, arity = 2, alphabet_size = 1):
		"""Creates a word object from a space separated string of letters. The *arity* and *alphabet_size* are automatically chosen to be as small as necessary, but they can be overridden if necessary. No :meth:`validation <validate>` is performed on the new word object.
		
			>>> w = Word.from_string("x a1 a1 a1 x a1 a1 a2 L a2")
			>>> print(w)
			x a1 a1 a1 x a1 a1 a2 L a2
		
		:raises ValueError: if the first symbol is not a letter :math:`x_i`.
		"""
		output = str.lower().split()
		if output[0][0] != 'x':
			raise ValueError("First symbol ({}) was not a letter x_i.".format(output[0]))
		for i, string in enumerate(output):
			char = string[0]
			if char in 'x':
				value = 1 if len(string) == 1 else int(string[1:])
				alphabet_size = max(alphabet_size, value)
				output[i] = value
			elif char in 'a':
				value = int(string[1:])
				arity = max(arity, value)
				output[i] = -value
			elif char in 'l': #ell for lambda
				output[i] = 0
		result = cls(arity, alphabet_size)
		result.extend(output)
		return result
	
	def __str__(self):
		"""Returns a string representation of this word. Passing the result into :meth:`from_string` will yield the original word. **NB:** words are not standardised before they are converted to a string."""
		return " ".join(self._char(i) for i in self)
	
	def _char(self, symbol):
		if symbol > 0:
			if self.alphabet_size == 1:
				return 'x'
			return 'x' + str(symbol)
		elif symbol < 0:
			return 'a' + str(-symbol)
		elif symbol == 0:
			return 'L'
	
	def __lt__(self, other):
		r"""Words are compared according to dictionary order (after reduction to :meth:`standard form <standardise>`). The alphabet :math:`X \cup \Omega` is ordered according to:
		
		..math :: \dots < \alpha_2 < \alpha_1 < \lambda < x_1 < x_2 < \dots
		"""
		#TODO examples and doctests
		if not isinstance(other, Word):
			return NotImplemented
		
		self.standardise()
		other.standardise()
		
		for s, o in zip(self, other):
			if s == o:
				continue
			if s < 0 and o < 0: #two alphas
				#self < o iff the index for s is less than the index for o
				return -s < -o
			return s < 0
		
		return False
	
	def __eq__(self, other):
		"""Words are standardised before testing for equality."""
		if not isinstance(other, Word):
			return NotImplemented
		self.standardise()
		other.standardise()
		return super().__eq__(self, other)
	
	def __ne__(self, other):
		if not isinstance(other, Word):
			return NotImplemented
		return not self.__eq__(other)
	
	def validate(self):
		"""Raises a ValueError if this word fails the valency test. (Prop 2.12)
		
			>>> w = Word.from_string("x a1 a2 L")
			>>> w.validate()
			Traceback (most recent call last):
				...
			ValueError: Word is invalid at symbol #4 (L).
			>>> w = Word.from_string("x a1 a2 x a2")
			>>> w.validate()
			Traceback (most recent call last):
				...
			ValueError: Word is invalid: valency is 2 (should be 1).
		"""
		valency = 0
		for i, symbol in enumerate(self):
			valency += self._valency_of(symbol)
			if valency <= 0:
				raise ValueError("Word is invalid at symbol #{} ({}).".format(i+1, self._char(symbol)))
		if valency != 1:
			raise ValueError("Word is invalid: valency is {} (should be 1).".format(valency))
	
	def _valency_of(self, symbol):
		"""Returns the valency of a symbol as per Def 2.11."""
		if symbol > 0: #letter x_i
			return 1
		if symbol < 0: #descenders alpha
			return 0
		if symbol == 0:#contractor lambda
			return 1 - self.arity
	
	@classmethod
	def initial_segment(cls, u, v):
		r"""Returns True if :math:`u` is an initial segment of :math:`v` or vice versa; otherwise returns False.
		
		Def 3.15: :math:`u` is an initial segment of :math:`v` if :math:`v = u\Gamma`, where :math:`\Gamma` is some string of alphas. Call :math:`u` a *proper* initial segment of :math:`v` if :math:`\Gamma` is not the empty string, *i.e.* if :math:`u \neq v`.
		
			>>> #completely different words
			>>> u = Word.from_string("x1 a2 a1 a1 a2")
			>>> v = Word.from_string("x2 a1 a2 a1 a2")
			>>> Word.initial_segment(u, v)
			False
			>>> #v starts with u
			>>> u = Word.from_string("x a1 a1 x2 a1 a2 L a1 a2")
			>>> v = Word.from_string("x a1 a1 a2 a2 a1")
			>>> Word.initial_segment(u, v)
			True
			>>> #v starts with u but has a lambda contraction later
			>>> u = Word.from_string("x a1 a2")
			>>> v = Word.from_string("x a1 a2 x a2 a2 L")
			>>> Word.initial_segment(u, v)
			False
		"""
		if len(u) > len(v):
			u, v = v, u
		u.standardise()
		v.standardise()
		return all(a == b for a, b in zip(u, v)) and all(b < 0 for b in v[len(u):])
		
			