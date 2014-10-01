"""
.. testsetup::
	
	from thompson.automorphism import *
"""

from collections import deque
from itertools import chain

from .word import Word, are_contractible
from .generators import Generators
from .full_tree import FullTree

def _concat(words):
	return list(chain.from_iterable(word for word in words)) + [0]

class Automorphism:
	r"""Represents an automorphism of :math:`V_{n,r}` by specifying two bases. This class keeps track of the mapping between bases.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	:ivar domain: a :class:`generating set <Generators>` of preimages 
	:ivar range: a :class:`generating set <Generators>` of images.
	"""
	
	def __init__(self, arity, alphabet_size, domain, range):
		r"""Creates an automorphism, given the *arity* :math:`n` and *alphabet_size* :math:`r`. Two bases *domain* and *range* are given. The automorphism maps elements so that order is preserved:
		
			.. math:: \text{domain}_i \mapsto \text{range}_i
		
		After creation, the automorphism is reduced. This means that collections of rules like
		
			.. math:: u\alpha_1 \mapsto v\alpha_1, \dotsc, u\alpha_n \mapsto v\alpha_n
		
		are replaced by simpler rules :math:`u \mapsto v` whenever possible.
		
		:raises ValueError: if the bases are of different sizes.
		:raises IndexError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis isn't actually a basis, i.e. is not a :meth:`free generating set <thompson.generators.Generators.is_free>` or does not :meth:`contract to the standard basis <thompson.generators.Generators.test_generates_algebra>`.
		"""
		
		#The boring checks
		if len(domain) != len(range):
			raise ValueError("Domain basis has {} elements, but range basis has {} elements.".format(
			  len(domain), len(range)))
		
		if not(arity == domain.arity == range.arity):
			raise IndexError("Arities do not match. Expected: {}, Domain: {}, Range: {}".format(
			  arity, domain.arity, range.arity))
		
		if not(alphabet_size == domain.alphabet_size == range.alphabet_size):
			raise IndexError("Alphabet sizes do not match. Expected: {}, Domain: {}, Range: {}".format(
			  alphabet_size, domain.alphabet_size, range.alphabet_size))
		
		#Check to see that the two generating sets given are free
		i, j = domain.test_free()
		if not(i == j == -1):
			raise ValueError("Domain is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		i, j = range.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		#Check to see that the two generating sets generate all of V_{n,r}
		missing = domain.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Domain does not generate V_{},{}. Missing elements are {}.".format(
			  arity, alphabet_size, [format(x) for x in missing]))
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{},{}. Missing elements are {}.".format(
			  arity, alphabet_size, [format(x) for x in missing]))
		
		#Before saving the domain and range, reduce them to remove any redundancy. This is like reducing tree pairs.
		Automorphism._reduce(domain, range)
		
		self.arity = arity
		self.alphabet_size = alphabet_size
		self.domain = domain
		self.range = range
		self._dict = {}
		for d, r in zip(self.domain, self.range):
			self._dict[d] = r
	
	@staticmethod
	def _reduce(domain, range):
		"""Contracts the domain generators whenever the corresponding contraction in range is possible. (This corresponds to reducing a tree pair diagram.)
			
			>>> from thompson.examples import cyclic_order_six as cyclic
			>>> for d, r in zip(cyclic.domain, cyclic.range):
			... 	print(d, '->', r)
			x1 a1 a1 -> x1 a1 a1
			x1 a1 a2 a1 -> x1 a1 a2 a2 a2
			x1 a1 a2 a2 -> x1 a2
			x1 a2 a1 -> x1 a1 a2 a2 a1
			x1 a2 a2 -> x1 a1 a2 a1
			>>> #Swaps x1 and x2
			>>> domain = Generators(2, 2, ["x1", "x2 a1", "x2 a2"])
			>>> range  = Generators(2, 2, ["x2", "x1 a1", "x1 a2"])
			>>> Automorphism._reduce(domain, range)
			>>> for d, r in zip(domain, range):
			... 	print(d, '->', r)
			x1 -> x2
			x2 -> x1
		"""
		#similar to word._reduce and Generator.test_generates_algebra
		i = 0
		arity = domain.arity
		while i <= len(domain) - arity:
			d_pref = are_contractible(domain[i : i + arity])
			r_pref = are_contractible(range[i : i + arity])
			if d_pref and r_pref: #are both non-empty tuples
				domain[i : i + arity] = [Word(d_pref, arity, domain.alphabet_size)]
				range[ i : i + arity] = [Word(r_pref, arity, range.alphabet_size )]
				i -= (arity - 1) 
				i = max(i, 0)
			else:
				i += 1
	
	def minimal_expansion(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \union W`. See example 4.25.
		
		>>> from thompson.examples import cyclic_order_six, example_4_25
		>>> cyclic_order_six.minimal_expansion()
		Generators(2, 1, ['x1 a1 a1', 'x1 a1 a2 a1', 'x1 a1 a2 a2', 'x1 a2'])
		>>> example_4_25.minimal_expansion()
		Generators(2, 1, ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
		"""
		basis = Generators.standard_basis(self.arity, self.alphabet_size)
		i = 0
		while i < len(basis):
			b = basis[i]
			if b in self.domain or b in self.range:
				i += 1
			else:
				basis.expand(i)
		return basis
	
	def to_quasinormal_form(self):
		"""The plan.
		
		1. Reduce the automorphism (eliminate carets) - DONE 
		2. Find elements above (Y union W) - DONE indirectly, see 3.
		3. Expand std basis until it contains all the elements from step 2. DONE
		4. Test each element of this expanded basis to see if the have nice orbits.
			a. If an element doesn't, expand it and try again with its children.
		5. When everything in the basis has nice orbits we are done.
		"""
		pass
	
	def __getitem__(self, word):
		"""Computes the image of a *word* under the given automorphism. The computation is cached for further usage later."""
		#Need to ensure this always returns a word
		#TODO:
			#if it's a tuple but not a word, validate
			#if it's a tuple and is a word, fine
			#if it's anything else, typeerror
		# if not isinstance(word, Word):
			# raise TypeError('{:r} is not a Word instance.'.format(word))
		#todo check word arity is the same and alphabet size not greater
		try:
			return self._dict[word]
		except KeyError:
			pass
		
		if word.is_simple():
			if self.domain.is_above(word):
				return self._image_of_simple(word)
			#otherwise, we have a simple word which is above our generators
			letters = _concat(self[child] for child in word.expand())
			image = Word(letters, self.arity, self.alphabet_size)
			self._dict[word] = image
			return image
		#else:
			# the word is complicated and in standard form, so ends with lambda
			# return Word(CONCAT(self[subword] for subword in lambda_arguments)
		return self.image_of_complex(word)
		
	def _image_of_simple(self, word):
		"""This method strips off alphas from the end of a string until the reamining string's image is already known. The computation is cached for further usage later.
		
			>>> from thompson.examples import example_4_25
			>>> #Basis element
			>>> u = Word('x1 a1', 2, 1)
			>>> print(example_4_25[u])
			x1 a1 a1 a1
			>>> #Below basis element
			>>> u = Word('x1 a1 a2 a2', 2, 1)
			>>> print(example_4_25[u])
			x1 a1 a1 a1 a2 a2
			>>> #Above basis element - no cancellation
			>>> u = Word('x1', 2, 1)
			>>> print(example_4_25[u])
			x1 a1 a1 a1 x1 a1 a1 a2 x1 a2 a2 x1 a1 a2 L x1 a2 a1 L L L
		"""
		i = 1
		while True:
			head, tail = word.split(i)
			if head in self._dict:
				break
			i += 1
			#should never loop forever here, provided the bases really are bases
			#and that everything above the bases has an image.
		
		image = self._dict[head]
		for _ in range(i):
			alpha, tail = tail[:1], tail[1:]
			head = head + alpha
			image = image.alpha(-alpha[0])
			self._dict[head] = image
		assert len(tail) == 0
		return image
		
	def image_of_complex(self, word):
		"""Computes the image of a non-simple *word*. This method forms a tree of lambda contractions, and works down the tree until it finds elements which we know the images of. The tree is then collapsed to form one contracted image.
		
		TODO doctest."""
		#Consume the lambda at the end of this word
		assert word[-1] == 0 #ends with a lambda
		root = FullTree(self.arity)
		root.data = [word, None]
		root.expand()
		for child in root:
			root.data = [GET_SUBWORDS_FROM_WORD, None]
		
		#data 0 and data 1 are clunky. change to pre image and post image, remove data from the tree class.
		
		unmapped = deque([root])
		while len(unmapped) > 0:
			node = unmapped.popleft()
			if node.data[0].is_simple():
				node.data[1] = self.image_of_simple(node.data[0]);
			else:
				node.expand()
				for child in node:
					child.data = [GET_SUBWORDS_FROM_WORD(node.data[0]), None]
					try:
						child.data[1]  = self._dict[child.data[0]]
					except KeyError:
						unmapped.append(child)
		
		#Now all nodes are mapped:
		for node in root.walk_postorder():
			if node.is_leaf():
				continue
			#TODO Word.concatenate
			node.data[1] = Word.concatenate((child.data[1] for child in node), self.arity, self.alphabet_size)
			self._dict[node.data[0]] = node.data[1]
		
		#Need to ensure that nodes above Y union W are put in the dict too. Use a similar contraction idea?
		
		
		
		
		
		