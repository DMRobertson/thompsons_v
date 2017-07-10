"""
.. testsetup::

	from thompson.word       import Word, format
	from thompson.generators import *
	from thompson.examples   import *
"""
from copy        import copy
from collections import defaultdict, Counter, Iterable

from .     import word
from .word import Word, Signature

__all__ = ["Generators", "rewrite_set", "rewrite_word"]

class Generators(list):
	r"""An ordered subset of :math:`V_{n, r}`, together with methods which treat such sets as generating sets and bases. Internally, this is a subclass of :class:`list <py3:list>`, so it is possible to reorder and otherwise modify collections of Generators.

	:ivar signature: The :class:`~thompson.word.Signature` of the algebra this set belongs to.
	"""

	def __init__(self, signature, generators=None):
		"""When creating a generating set, you must specify the algebra :math:`V_{n, r}` it is a subset of. A list of generators can optionally be provided. Each generator is passed to :meth:`append`.
		"""
		if not isinstance(signature, Signature):
			signature = Signature(*signature)
		self.signature = signature

		self.cache = None

		if generators is not None:
			for g in generators:
				self.append(g)

	def append(self, w):
		"""Adds *w* to this generating set. If *w* is a string, a :class:`~thompson.word.Word` object is created and assigned the same *arity* and *alphabet_size* as the generating set.

		:raises TypeError: if *w* is neither a string nor a Word.
		:raises IndexError: if *w* has a different arity to the generating set.
		:raises IndexError: if *w* has a larger alphabet_size the generating set.
		:raises ValueError: if *w* is already contained in this generating set.
		"""
		if isinstance(w, str):
			w = Word(w, self.signature)
		elif not isinstance(w, Word):
			raise TypeError("{!r} is neither a string nor a Word.".format(w))

		if self.signature != w.signature:
			raise IndexError("Can't add {} with arity {} to generating set with arity {}."
			  .format(w, w.signature, self.signature))

		if w not in self:
			super().append(w)
		else:
			raise ValueError("The word {} is already in this set.".format(w))
	
	def extend(self, words):
		#This is probably inefficient but I just want an extend method that does the checks that append() does above
		for word in words:
			self.append(word)
	
	def __add__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		new = self.copy()
		new.extend(other)
		return new
		
	def __iadd__(self, other):
		if not isinstance(other, Iterable):
			return NotImplemented
		self.extend(other)

	def __str__(self):
		return "[" + ", ".join(format(w) for w in self) + "]"
	
	def __format__(self, specification):
		return str(self).__format__(specification)
	
	def __repr__(self):
		return type(self).__name__ + "({}, [{}])".format(
		  self.signature, ", ".join(repr(format(w)) for w in self))
	
	def __getitem__(self, key):
		if isinstance(key, slice):
			output = type(self)(self.signature)
			output.extend( super().__getitem__(key) )
			return output
		return super().__getitem__(key)
	
	def __eq__(self, other):
		"""Two Generators instances are equal iff they have the same signature and the same elements in the same order.

			>>> x = random_basis(signature=(3, 2))
			>>> y = x.embed_in((4, 3))
			>>> x == y #Different signatures, so False
			False
			>>> list.__eq__(x, y) #Even though they have the same elements
			True
		"""
		if not isinstance(other, Generators):
			return NotImplemented
		return self.signature == other.signature and super().__eq__(other)

	def copy(self):
		"""We override ``list.copy()`` so that we don't have to recreate :class:`Generators` instances by hand.

			>>> olga_f = load_example('olga_f')
			>>> X = olga_f.quasinormal_basis.copy()
			>>> X is olga_f.quasinormal_basis
			False
			>>> type(X).__name__
			'Generators'
		"""
		return type(self)(self.signature, self)

	def filter(self, func):
		"""Creates a copy of the current generating set whose elements *x* are those for which *func(x)* is ``True``. The original generating set is unmodified, and the order of elements is inherited by the filtered copy.

			>>> X = load_example('olga_f').quasinormal_basis
			>>> print(X)
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2]
			>>> def test(x):
			... 	return len(x) % 2 == 0
			...
			>>> print(X.filter(test))
			[x1 a2 a2 a1, x1 a2 a2 a2]
		"""
		return type(self)(self.signature, (x for x in self if func(x)))

	#Tests on generating sets
	def test_free(self):
		r"""Tests to see if the current generating set :math:`X` is a free generating set.
		To do so, the words contained in :math:`X` must all be :meth:`simple <thompson.word.Word.is_simple>`.

		If the test fails, returns the first pair of indices :math:`(i, j)` found for which :math:`X_i` :meth:`is above <thompson.word.Word.is_above>` :math:`X_j`. If the test passes, returns :math:`(-1, -1)`.

			>>> g = Generators((2, 3), ["x1 a1", "x1 a2 a1", "x1 a2 a1 a1", "x1 a2 a2"])
			>>> g.test_free()
			(1, 2)
			>>> print(g[1], g[2], sep='\n')
			x1 a2 a1
			x1 a2 a1 a1
			>>> g = Generators((2, 1), ['x1 a2 a2', 'x1 a2 a1 x1 a2 a2 a2 a1 L x1 a1 L'])
			>>> g.test_free()
			Traceback (most recent call last):
			...
			ValueError: Cannot test for freeness unless all words are simple.

		:raises ValueError: if any word in the basis is not simple.

		.. seealso::  Lemma :paperref:`HigmanLemma2.5` of the paper.
		"""
		if any(not word.is_simple() for word in self):
			raise ValueError('Cannot test for freeness unless all words are simple.')

		for i in range(len(self)):
			for j in range(i+1, len(self)):
				if self[i].is_above(self[j]):
					return i, j
				elif self[j].is_above(self[i]):
					return j, i
		return -1, -1

	def is_free(self):
		"""Returns True if this is a free generating set; otherwise False.

			>>> g = Generators((2, 3), ['x1', 'x2']);
			>>> g.is_free()
			True
			>>> g.append('x3'); g.is_free()
			True
			>>> g.append('x2 a1'); g.is_free()
			False
		"""
		return self.test_free() == (-1, -1)

	def preorder_traversal(self):
		"""Yields words as if traversing the vertices of a tree in pre-order."""
		if self.signature.alphabet_size > 1:
			raise NotImplementedError('Not implemented for alphabet sizes greater than 1 ({})'.format(
			  self.signature.alphabet_size))

		yield Word('x', self.signature), (None if len(self) > 1 else 1)
		path = Word('x a1', (2, 1))
		while len(path) > 1:
			leaf = path in self
			leafnum = self.index(path) + 1 if leaf else None
			yield path, leafnum
			if leaf:
				finished = False
				while not finished:
					last = abs(path[-1])
					path = Word(path[:-1], self.signature)
					finished = last != self.signature.arity or len(path) == 1
				if last != self.signature.arity:
					path = path.alpha(last + 1)
			else:
				path = path.alpha(1)

	def simple_words_above(self):
		r"""An iterator that yields the simple words which can be obtained by contracting this basis :math:`n \geq 0` times. (So the 'above' condition is not strict.)

			>>> g = Generators((3, 1), ["x a1 a1 a3", "x a1 a1 a1", "x a1 a2", "x a1 a1 a2", "x a1 a3"])
			>>> for word in g.simple_words_above():
			... 	print(format(word))
			x1 a1 a1 a1
			x1 a1 a1 a2
			x1 a1 a1 a3
			x1 a1 a2
			x1 a1 a3
			x1 a1 a1
			x1 a1
		"""
		words = sorted(self) #Creates a shallow copy and then sorts it.
		yield from words
		i = 0
		arity = self.signature.arity
		while i <= len(words) - arity:
			prefix = word.are_contractible(words[i : i + arity])
			if prefix:
				yield prefix
				words[i : i + arity] = [prefix] if prefix not in words else []
				i -= (arity - 1)
				i = max(i, 0)
			else:
				i += 1

	def test_generates_algebra(self):
		r"""Tests to see if this set generates all of :math:`V_{n,r}`. The generating set is sorted, and then contracted as much as possible. Then we check to see which elements of the :meth:`standard basis <standard_basis>` :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` are not included in the contraction.

		The test fails if there is at least one such element; in this case, the list of all such elements is returned. Otherwise the test passes, and an empty list is returned.

			>>> #A basis for V_2,1
			>>> Y = Generators((2, 1), ["x a1", "x a2 a1", "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
			>>> Y.test_generates_algebra()
			[]
			>>> #The same words do not generate V_2,3
			>>> Y = Generators((2, 3), ["x a1", "x a2 a1", "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
			>>> missing = Y.test_generates_algebra(); missing
			[Word('x2', (2, 3)), Word('x3', (2, 3))]
			>>> for x in missing: print(x)
			x2
			x3

		:rtype: a list of :class:`Words <thompson.word.Word>`."""
		words = set(self.simple_words_above())
		#At the end, should contract to [x1, x2, ..., x_r]
		expected = Generators.standard_basis(self.signature)
		missing = [x for x in expected if x not in words]
		return missing

	def generates_algebra(self):
		"""Returns True if this set generates all of :math:`V_{n,r}`. Otherwise returns False.

			>>> from random import randint
			>>> arity = randint(2, 5); alphabet_size = randint(2, 10)
			>>> Y = Generators.standard_basis((arity, alphabet_size))
			>>> Y.generates_algebra()
			True
			>>> random_basis().generates_algebra()
			True

		"""
		return len(self.test_generates_algebra()) == 0

	def test_above(self, word, return_index=False):
		"""Searches for a pair (gen, tail) where gen belongs to the current basis and gen + tail = word. If no such pair exists, returns None.

		If *return_index* is False, returns the pair *(gen, tail)*. Otherwise, returns the pair *(i, tail)* where *i* is the index of *gen* in the current basis.

			>>> basis = Generators.standard_basis((2, 1)).expand(0).expand(0).expand(0)
			>>> basis
			Generators((2, 1), ['x1 a1 a1 a1', 'x1 a1 a1 a2', 'x1 a1 a2', 'x1 a2'])
			>>> gen, tail = basis.test_above(Word('x1 a2 a2 a1', (2, 1)))
			>>> print(gen, '|', format(tail))
			x1 a2 | a2 a1
			>>> i, tail = basis.test_above(Word('x1 a2 a2 a1', (2, 1)), return_index=True)
			>>> print(basis[i])
			x1 a2
			>>> basis.test_above(Word('x1', (2,1))) is None
			True
			>>> gen, tail = basis.test_above(basis[0])
			>>> print(gen)
			x1 a1 a1 a1
			>>> gen is basis[0] and len(tail) == 0
			True
		"""
		if not return_index and self.cache is not None:
			return self._test_above_cached(word)
		for i, gen in enumerate(self):
			result = gen.test_above(word)
			if result is not None:
				if return_index:
					return i, result
				return gen, result
		return None

	def _test_above_cached(self, word):
		"""A quicker version of :meth:`test_above` which assumes that the:
			- a cache set has been maintained for this generating set, and
			- each generator is simple.

		.. todo::

			The bulk of the speedup seemed to be properly using all the information available. Do we still need this cache? Given the fact that the QNB is used virtually all the time it might be a good idea.
		"""
		if not word.is_simple():
			return None
		for i in range(len(word)+1):
			if word[:i] in self.cache:
				head, tail = word.split(i)
				return Word(head, word.signature, preprocess=False), tail
		return None

	def is_above(self, word):
		"""Returns True if any generator :meth:`~thompson.word.Word.is_above` the given *word*.

			>>> g = Generators((2, 2), ['x1 a1', 'x1 a2', 'x2'])
			>>> g.is_above(Word('x1 a1 a1 a2', (2, 2)))
			True
			>>> g.is_above(Word('x1', (2, 2)))
			False
			>>> g.is_above(Word('x2', (2, 2)))
			True
			>>> g.is_above(Word('x1 a1 x1 a2 L', (2, 2)))
			False
		"""
		return any(gen.is_above(word) for gen in self)

	def is_basis(self):
		"""Returns True if this is a :meth:`free generating set <is_free>` which :meth:`generates all of <is_free>` :math:`V_{n,r}`. Otherwise returns False.

			>>> g = Generators((2, 2), ["x1 a1", "x1 a2"])
			>>> g.is_free()
			True
			>>> g.generates_algebra()
			False
			>>> g.is_basis()
			False
			>>> g.append('x2')
			>>> g.is_basis()
			True
			>>> random_basis().is_basis()
			True
		"""
		return self.is_free() and self.generates_algebra()

	#Information about generating sets
	def descendants_above(self, floor):
		r"""Suppose we have a basis *floor* below the current basis *ceiling*. There are a finite number of elements below *ceiling* which are not below *floor*. This method enumerates them. In symbols, we are enumerating the set :math:`X\langle A\rangle \setminus F\langle A\rangle`, where :math:`X` is the current basis and :math:`F` is the *floor*.

			>>> phi = load_example('example_5_15')
			>>> X = phi.quasinormal_basis
			>>> Y = X.minimal_expansion_for(phi)
			>>> Z = phi.image_of_set(Y)
			>>> terminal = X.descendants_above(Y)
			>>> initial  = X.descendants_above(Z)
			>>> print(initial, terminal, sep='\n')
			[x1 a1, x1 a1 a1]
			[x1 a2 a2, x1 a2 a2 a1]
		"""
		if not all(self.is_above(f) for f in floor):
			raise ValueError('Floor is not below the current generating set.')

		output = copy(self)
		i = len(output) - 1
		while i >= 0:
			if output[i] in floor:
				del output[i]
			i -= 1

		i = 0
		while i < len(output):
			word = output[i]
			children = [child for child in word.expand() if child not in floor]
			output.extend(children)
			i += 1
		return output

	#Creating generating sets
	@classmethod
	def standard_basis(cls, signature):
		r"""Creates the standard basis :math:`\boldsymbol{x} = \{x_1, \dotsc, x_r\}` for :math:`V_{n,r}`, where :math:`n` is the arity and :math:`r` is the *alphabet_size* of the *signature*.

			>>> Generators.standard_basis((2, 4))
			Generators((2, 4), ['x1', 'x2', 'x3', 'x4'])
		"""
		if isinstance(signature, tuple) and not isinstance(signature, Signature):
			signature = Signature(*signature)

		output = cls(signature)
		for i in range(1, signature.alphabet_size + 1):
			output.append(Word((i,), signature, preprocess=False))
		return output
	
	@classmethod
	def from_dfs(cls, string):
		"""Creates bases correponding to binary trees via a string of 1s and 0s correponding to branches and carets respectively. See :meth:`~thompson.automorphism.Automorphism.from_dfs`. The *string* can also be specified as a :class:`~py3:int`, in which case it is replaced by its base 10 string representation.
		
			>>> print(Generators.from_dfs("100"))
			[x1 a1, x1 a2]
			>>> print(Generators.from_dfs(1100100)) #ints are fine too
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
		
		:raises ValueError: if the *string* doesn't correctly describe a rooted binary tree.
		
		.. doctest::
			
			>>> print(Generators.from_dfs(""))
			Traceback (most recent call last):
			...
			ValueError: Should have one more zero than one (got 0 and 0 respectively)
			>>> print(Generators.from_dfs("10"))
			Traceback (most recent call last):
			...
			ValueError: Should have one more zero than one (got 1 and 1 respectively)
			>>> print(Generators.from_dfs(10001))
			Traceback (most recent call last):
			...
			ValueError: Error at character 3: complete description of tree with unprocessed digits remaining. String was 10001
		"""
		if isinstance(string, int):
			string = str(string)
		counts = Counter(string)
		if not set(counts) <= set('01 '):
			raise ValueError("String {} should contain only 1s, 0s and optionally spaces".format(string))
		if counts['0'] != counts['1'] + 1:
			raise ValueError("Should have one more zero than one (got {} and {} respectively)".format(
				counts['0'], counts['1']))
		
		basis = cls.standard_basis((2, 1))
		index = 0
		counts = Counter()
		for string_index, char in enumerate(string, start=1):
			if char == ' ':
				continue
			counts[char] += 1
			balance = counts['0'] - (counts['1'] + 1)
			if balance > 0:
				raise ValueError("Error at character {}: too many zeroes before a one. String was {}".format(
					string_index, string))
			elif balance == 0 and string_index != len(string):
				#should now be at the end of the string. If not, error
				raise ValueError("Error at character {}: complete description of tree with unprocessed digits remaining. String was {}".format(
					string_index, string))
			
			#else all is well, handle this char.
			if char == '1':
				basis.expand(index)
			elif char == '0':
				index += 1
		
		return basis
	
	#TODO: Generators.to_dfs(self):
		#...
	
	def minimal_expansion_for(self, *automorphisms):
		r"""Suppose we are given a finite sequence of *automorphisms* of :math:`V_{n,r}` and that the current generaeting set is a basis :math:`X` for :math:`V_{n,r}`. This methods returns an expansion :math:`Y` of :math:`X` such that each automorphism maps :math:`Y` into :math:`X\langle A\rangle`.

			>>> std_basis = Generators.standard_basis((2, 1))
			>>> reduced_domain = Generators((2, 1), ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2", "x a2 a1", "x a2 a2"])
			>>> std_basis.minimal_expansion_for(load_example('cyclic_order_six')) == reduced_domain
			True

			>>> phi = load_example('example_5_3')
			>>> basis = phi.quasinormal_basis
			>>> print(basis.minimal_expansion_for(phi))
			[x1 a1 a1 a1 a1, x1 a1 a1 a1 a2, x1 a1 a1 a2, x1 a1 a2 a1, x1 a1 a2 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2]

		:raises ValueError: if no automorphisms are passed to this method.
		:raises ValueError: if the automorphisms or basis do not all belong to the same group :math:`G_{n,r}`.
		:raises ValueError: if the given basis does not generate :math:`V_{n,r}`.

		.. seealso:: Lemma :paperref:`4.1H` of the paper proves that this expansion exists, is of minimal size, and is unique with this that property.
		"""
		#0. Pre-flight checks.
		if len(automorphisms) == 0:
			raise ValueError('Must provide at least one automorphism to minimal_expansion_for().')

		basis = copy(self)
		if not basis.is_free():
			raise ValueError('The generating set does not freely generate V_{n,r}.')

		for aut in automorphisms:
			if aut.signature != self.signature:
				raise ValueError('At least one automorphism belongs to a different G_{n,r} than the basis.')

		#1. Expand until all images belong to X<A>.
		i = 0
		while i < len(basis):
			b = basis[i]
			if all(self.is_above(aut.image(b)) for aut in automorphisms):
				i += 1
			else:
				basis.expand(i)
		return basis

	def embed_in(self, signature, shift=0):
		"""Creates a copy of the current generating set in the algebra determined by *signature*.

		:raises ValueError: if the current signature's algebra is not a subset of the given signature's algebra.

			>>> x = Generators.standard_basis((2, 1)).expand(0); x
			Generators((2, 1), ['x1 a1', 'x1 a2'])
			>>> y = x.embed_in((3, 2)); y
			Generators((3, 2), ['x1 a1', 'x1 a2'])
			>>> z = x.embed_in((3, 2), shift=1); z
			Generators((3, 2), ['x2 a1', 'x2 a2'])
			>>> print(y.signature, y[0].signature)
			(3, 2) (3, 2)
			>>> x.is_basis()
			True
			>>> y.is_basis()
			False

		"""
		if not isinstance(signature, Signature):
			signature = Signature(*signature)
		if (self.signature.arity > signature.arity or
		    self.signature.alphabet_size + shift > signature.alphabet_size):
			raise ValueError("Cannot embed {} and shift by {} into {}".format(
			  self.signature, shift, signature))

		output = Generators(signature)
		for word in self:
			output.append( word.shift(shift, signature) )
		return output

	#Modifiers
	def expand(self, index):
		r"""Replaces the word :math:`w` at index *index* in the current generating set with its children, :math:`w\alpha1, \dotsc, w\alpha_n`, where :math:`n` is the arity of the generating set. As with ordinary Python lists, negative values of *index* are supported.

			>>> g = Generators.standard_basis((3, 1)); g
			Generators((3, 1), ['x1'])
			>>> g.expand(0)
			Generators((3, 1), ['x1 a1', 'x1 a2', 'x1 a3'])
			>>> g #has been modified
			Generators((3, 1), ['x1 a1', 'x1 a2', 'x1 a3'])
			>>> g.expand(-2) #expand at the second entry from the right, i.e. at 'x1 a2'
			Generators((3, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2', 'x1 a2 a3', 'x1 a3'])

			:raises IndexError: if there is no generator at index *index*.
			:returns: the current generating set, after modification.
		"""
		if -len(self) <= index < 0:
			index %= len(self)
		elif index < 0:
			raise IndexError('Index {} out of range 0..{} or -{}..-1'.format(index, len(self), -len(self)))
		self[index: index+1] = [self[index].alpha(i) for i in range(1, self.signature.arity+1)]
		return self #allows chaining

	def _expand_with_cache(self, index):
		r"""If a cache is being maintained, this function expands() and maintains the cache."""
		assert self.cache is not None
		self.cache.remove(self[index])
		self.expand(index)
		self.cache.update(self[index : index + self.signature.arity])
		return self
	
	def expand_away_lambdas(self):
		r"""Experimental method to recursively (and inefficiently) expand a generating set until every word has lambda-length 0.
		
			>>> X = Generators((2,1), ['x1 a1 a1 a1 a2', 'x1 a2 a1 a2 x1 a2 a1 a1 a2 x1 a2 a2 a2 L x1 a2 a2 a1 L L'])
			>>> X.expand_away_lambdas()
			>>> print(X)
			[x1 a1 a1 a1 a2, x1 a2 a1 a2, x1 a2 a1 a1 a2, x1 a2 a2 a2, x1 a2 a2 a1]
		"""
		index = 0
		while index < len(self):
			if self[index].lambda_length > 0:
				self.expand(index)
			else:
				index += 1
		
	@classmethod
	def sort_mapping_pair(cls, domain, range):
		"""Makes copies of the given lists of words *domain* and *range*. The copy of *domain* is sorted according to the :meth:`order <thompson.word.Word.__lt__>` defined on words. The same re-ordering is applied to the *range*, so that the mapping from *domain* to *range* is preserved.

		:rtype: A pair of class:`Generator` instances.
		"""
		d, r = zip(*sorted(zip(domain, range)))
		return cls(domain.signature, d), cls(range.signature, r)

	def expand_to_size(self, size):
		"""Expands the current generating set until it has the given *size*. The expansions begin from the end of the generating set and work leftwards, wrapping around if we reach the start. (This is to try and avoid creating long words where possible.)
		
			>>> basis = Generators.standard_basis((3, 1)); print(basis)
			[x1]
			>>> basis.expand_to_size(11); print(basis)
			[x1 a1 a1, x1 a1 a2, x1 a1 a3, x1 a2 a1, x1 a2 a2, x1 a2 a3, x1 a3 a1, x1 a3 a2, x1 a3 a3 a1, x1 a3 a3 a2, x1 a3 a3 a3]
			>>> basis = Generators.standard_basis((2, 1)); print(basis)
			[x1]
			>>> basis.expand_to_size(2); print(basis)
			[x1 a1, x1 a2]
			>>> basis = Generators.standard_basis((2, 3)).expand(2).expand(0); print(basis)
			[x1 a1, x1 a2, x2, x3 a1, x3 a2]
			>>> basis.expand_to_size(12); print(basis)
			[x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2, x2 a1, x2 a2, x3 a1 a1, x3 a1 a2, x3 a2 a1 a1, x3 a2 a1 a2, x3 a2 a2 a1, x3 a2 a2 a2]
			
		Expanding a basis to its current size does nothing.
		
			>>> b1 = random_basis()
			>>> b2 = b1.copy()
			>>> b2.expand_to_size(len(b1))
			>>> b1 == b2
			True
		
		If expansion to the target size is not possible, a :class:`py3:ValueError` is raised.
		
			>>> basis = Generators.standard_basis((3,2))
			>>> len(basis)
			2
			>>> basis.expand_to_size(3)
			Traceback (most recent call last):
			...
			ValueError: Cannot expand from length 2 to length 3 in steps of size 2.
			>>> basis.expand_to_size(4)
			>>> len(basis)
			4
			>>> basis.expand_to_size(1)
			Traceback (most recent call last):
			...
			ValueError: Cannot expand from length 4 to length 2 in steps of size 2.
		"""
		modulus = self.signature.arity - 1
		if (size % modulus != len(self) % modulus) or size < len(self):
			raise ValueError("Cannot expand from length {} to length {} in steps of size {}.".format(
			  len(self), size, modulus))

		num_expansions = (size - len(self))// modulus
		i = -1
		for _ in range(num_expansions):
			if i < 0:
				i = len(self) - 1
			self.expand(i)
			i -= 1
		assert len(self) == size

	def cycled(self, positions=1):
		"""Produces a copy of the current generating set, cyclicly shifted by a number of *positions* to the right.

			>>> phi = load_example('example_5_15')
			>>> print(phi.domain)
			[x1 a1, x1 a2 a1, x1 a2 a2 a1 a1, x1 a2 a2 a1 a2, x1 a2 a2 a2]
			>>> print(phi.domain.cycled(2))
			[x1 a2 a2 a1 a1, x1 a2 a2 a1 a2, x1 a2 a2 a2, x1 a1, x1 a2 a1]

			>>> X = random_basis()
			>>> X == X.cycled(0) == X.cycled(len(X))
			True
		"""
		positions = positions % len(self)
		out = type(self)(self.signature)
		for i in range(positions, len(self)):
			out.append(self[i])
		for i in range(positions):
			out.append(self[i])
		return out

def rewrite_set(set, basis, new_basis):
	"""Maps *set* into the algebra generated by *new_basis* according to the bijection between *basis* and *new_basis*. See also :func:`rewrite_word`."""
	set.signature = new_basis.signature
	for i, word in enumerate(set):
		set[i] = rewrite_word(word, basis, new_basis)

def rewrite_word(word, basis, new_basis):
	r"""Suppose that we have a *word* which is below a given *basis*. Suppose we have a bijection between *basis* and some image *new_basis*. Then we can rewrite *word* as as descendant of *new_basis* in a manner which is compatible with this bijection.

	:raises ValueError: if *basis* is not above *word*.
	:raises ValueError: if *basis* and *new_basis* have different sizes.
	"""
	if len(basis) != len(new_basis):
		raise ValueError('Basis sizes {} and {} do not match.'.format(
		  len(basis), len(new_basis)))

	result = basis.test_above(word, return_index=True)
	if result is None:
		raise ValueError('The word {} is not below the basis {}'.format(
		  word, basis))
	index, tail = result
	return new_basis[index].extend(tail)
