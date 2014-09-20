"""
.. testsetup:: 
	
	from thompson.permutation import *
"""

from fractions import gcd
from functools import reduce

__all__ = ["Permutation"]

class Permutation:
	r"""A permutation *f* of a finite set :math:`\{1, \dotsc, N\}` of integers.
	
	:ivar size: the size :math:`N` of the set being permuted.
	"""
	def __init__(self, output):
		r"""This method accepts the *output* :math:`[f(1), f(2), \dotsc, f(N)]` of *f*, specified as
		
		- a list of integers, or
		- a space-separated string of base 10 integers (which is converted into the above).
		
		Put differently, the input is the bottom row of a `two-row matrix representation <http://en.wikipedia.org/wiki/Permutation#Definition_and_usage>`_
		
		.. math::
			
			\begin{pmatrix}
				1		&2		&\dots	&N	\\
				f(1)	&f(2)	&\dots	&f(N)
			\end{pmatrix}
		
		of *f*.
		
		:raises ValueError: if the given output does not properly define a permutation. For example:
		
			>>> Permutation('1 2 1 3') #reused 1
			Traceback (most recent call last):
			...
			ValueError: Label 1 reused in [1, 2, 1, 3]
			>>> Permutation('1 5 2 7') #doesn't map {1, 2, 3, 4} to itself
			Traceback (most recent call last):
			...
			ValueError: Elements [3, 4] are not included in the output [1, 5, 2, 7]
		"""
		if isinstance(output, str):
			temp = [int(x) for x in output.split()]
			assert len(temp) > 0, "Permutation created from {} with of length zero".format(repr(output))
			output = temp
		
		#1. Check that the output don't contain any duplicates.
		alphabet = set()
		for item in output:
			if item in alphabet:
				raise ValueError("Label %r reused in %r" % (item, output))
			alphabet.add(item)
		size = len(output)
		
		#2. Check that this really is a bijection of {1, ...., size}
		domain = { i + 1 for i in range(size)}
		missing = domain - alphabet
		if len(missing) > 0:
			raise ValueError("Elements %s are not included in the output %r."
			  % (sorted(missing), output ) )
		
		#3. Set up the useful stuff
		self.output = output
		self.size = size
	
	def __getitem__(self, key):
		#TODO. Might be more natural for a mathematician to use __call__? perm(1) instead of perm[1]?
		r"""Permutations can be accessed like a dictionary or list to map inputs to output. The image :math:`f(x)` of  :math:`x \in \{1, \dotsc, N\}` is given by ``f[x]``.
		
			>>> f = Permutation("4 3 1 2")
			>>> f[1]
			4
			>>> f[10]
			Traceback (most recent call last):
			...
			KeyError: "Key 10 is out of range 1, ..., 4."
		"""
		if not (1 <= key <= self.size):
			raise KeyError("Key %i is out of range 1, ..., %i." % (key, self.size))
		return self.output[key-1]
	
	def __setitem__(self, key, value):
		"""Modifies the permutation so that the image of *key* under the permutation is *value*. Care should be taken to ensure that the permutation is still a permutation after modification.
		
		:raises ValueError: if *value* is not in the range 1...n, where n is the permutation's size.
		:raises KeyError: if the same is true of *key*.
		
			>>> x = Permutation("3 2 1 4") #(1 3)
			>>> x[10] = 8
			Traceback (most recent call last):
				...
			KeyError: Key 10 is out of range 1, ..., 4.
			>>> x[1] = 8
			Traceback (most recent call last):
				...
			ValueError: Value 8 is out of range 1, ..., 4.
			>>> x[1] = 1
			>>> x[3] = 3
			>>> print(x) #x is now the identity
			()
		"""
		if not (1 <= key <= self.size):
			raise KeyError("Key %i is out of range 1, ..., %i." % (key, self.size))
		if not (1 <= value <= self.size):
			raise ValueError("Value %i is out of range 1, ..., %i." % (value, self.size))
		self.output[key-1] = value
	
	def __iter__(self):
		"""Iterating over a permutation yields tuples (x, y) where y is the image of x under the current permutation.
		
			>>> for x, y in Permutation("2 3 1 4"):
			... 	print(x, '->', y)
			1 -> 2
			2 -> 3
			3 -> 1
			4 -> 4
		"""
		for i in range(self.size):
			yield (i+1), self.output[i]
	
	def __str__(self):
		"""When permutations are converted to a string, they are represented as a product of cycles.
		
			>>> str(Permutation("2 3 1 5 4"))
			'(1 2 3)(4 5)'
		"""
		cycles = self.cycles()
		if len(cycles) == 0:
			return "()"
		out = ""
		for cycle in cycles:
			out += "(" + " ".join(str(x) for x in cycle) + ")"
		return out
	
	def __repr__(self):
		out_str = " ".join(str(x) for x in self.output)
		return "%s(%r)" % (self.__class__.__name__, out_str)
	
	def __eq__(self, other):
		"""Two permutations are equal if they have the same size and 
		
			>>> Permutation("2 3 1") == Permutation("3 2 1")
			False
			>>> Permutation("1 2 3") == Permutation("1 2 3 4") #Identities of S_3 and S_4
			False
			>>> x = Permutation("2 1 4 3") #Product of 2-cycles (1 2)(3 4)
			>>> x == x.inverse()
			True
		"""
		if not isinstance(other, Permutation):
			return NotImplemented
		return self.output == other.output
	
	def __ne__(self, other):
		if not isinstance(other, Permutation):
			return NotImplemented
		return not(self.output == other.output)
	
	def is_identity(self):
		"""Returns ``True`` if the permutation maps *x* to itself, for all applicable *x*. Otherwise ``False``.
		
			>>> Permutation("1 2 3 4").is_identity()
			True
			>>> Permutation("2 1 3 4").is_identity()
			False
		"""
		return all(x == image for (x, image) in self)
	
	def is_cycle(self, of_length=None):
		"""Returns ``True`` if the permutation consists of a single cycle, otherwise ``False``. 
		
			>>> Permutation("1 2 3 4").is_cycle() #identity
			True
			>>> Permutation("2 1 3 4").is_cycle() #2-cycle (1 2)
			True
			>>> Permutation("2 1 4 3").is_cycle() #product of 2-cycles (1 2)(3 4)
			False
		
		If the argument *of_length* is given, the method returns ``True`` only if the permutation consists of a single cycle whose length is *of_length*. Identity permutations count as having length 0 or 1 (though to test for identity permutations it's better to use :meth:`is_identity`.
		
			>>> Permutation("1 2 3 4").is_cycle(of_length=0) #identity
			True
			>>> Permutation("1 2 3 4").is_cycle(of_length=1) #identity
			True
			>>> Permutation("2 1 3 4").is_cycle(of_length=2) #2-cycle (1 2)
			True
			>>> Permutation("2 1 3 4").is_cycle(of_length=3) #2-cycle (1 2)
			False
			>>> Permutation("2 1 4 3").is_cycle(of_length=2) #product of 2-cycles (1 2)(3 4)
			False
		"""
		cycles = self.cycles()
		num_cycles = len(cycles)
		#1. Return false if this is a product of cycles.
		if num_cycles > 1:
			return False
		#2. Identity permutations are single cycles; check to if this hasthe required length.
		if num_cycles == 0:
			return of_length in (0, 1, None)
		#3. Now we have just a single cycle which is not the identity.
		return (of_length is None) or (len(cycles[0]) == of_length)
	
	def cycles(self):
		r"""Returns a list of tuples representing this permutation as a product of cycles. Trivial cycles such as (6), (3) and () correspond to the identity permutation, and are not included in the output.
		
			>>> Permutation("2 3 4 5 6 1").cycles()
			[(1, 2, 3, 4, 5, 6)]
			>>> Permutation("2 3 1 5 4 6").cycles()
			[(1, 2, 3), (4, 5)]
		"""
		out = list(range(1, self.size + 1))
		cycles = []
		while out:
			start = out.pop(0)
			orbit = [start]
			next = self[start]
			while next != start:
				orbit.append(next)
				out.remove(next)
				next = self[next]
			if len(orbit) > 1:
				cycles.append(tuple(orbit))
		
		return cycles
	
	def order(self):
		"""Returns the order of the permutation. This is the smallest positive integer *n* such that the *n* th power of the permutation is the identity. Identity permutations have an order of 1.
		
			>>> Permutation("1 2 3 4 5 6").order() #Identity
			1
			>>> Permutation("2 3 4 5 6 1").order() #6-cycle (1 2 3 4 5 6)
			6
			>>> Permutation("2 1 3 4 6 5").order() #product of 2-cycles (1 2)(5 6)
			2
			>>> Permutation("2 3 1 4 6 5").order() #product of cycle (1 2 3)(5 6)
			6
		"""
		return lcm(len(cycle) for cycle in self.cycles())
	
	def __mul__(self, other): #self after other
		"""Multiplying permutations yields their composition.
		
			>>> f = Permutation("1 2 4 3") #(3 4)
			>>> g = Permutation("3 2 1 4") #(1 3)
			>>> print(f * g) #f after g
			(1 4 3)
			>>> a = Permutation("2 1") #(2 1)
			>>> b = Permutation("5 1 2 3 4") #(1 5 4 3 2)
			>>> print(a * b) #a after b
			(1 5 4 3)
			>>> print(b * a) #b after a
			(2 5 4 3)
		"""
		size = max(self.size, other.size)
		output = [0] * size
		for i in range(size):
			try:
				x = other[i+1]
			except KeyError:
				x = i+1
			try:
				x = self[x]
			except KeyError:
				x = x
			output[i] = x
		return type(self)(output)
	
	def inverse(self):
		"""Creates and returns a new permutation object representing the inverse of the current permutation.
		
		>>> from random import randint, shuffle
		>>> x = list(range(1, randint(5,10)))
		>>> shuffle(x)
		>>> p = Permutation(x)
		>>> print(p)
		(...)
		>>> q = p.inverse()
		>>> print(q)
		(...)
		>>> (p * q).is_identity()
		True
		"""
		output = [0] * self.size
		for i, image in self:
			output[image-1] = i
		return type(self)(output)
	
	def remove_image_of(self, index):
		r"""Removes the rule that index -> self[index] from the permutation. The permutation is relabelled so that it uses the symbols 1 to n-1 instead of 1 to n.
		
		**Example**. Suppose we begin with ``Permutation("3 4 5 6 1 2")`` and suppose we wish to remove the image of 2. We delete the rule 2 -> 4 (coloured red). Next, we reduce the domain labels (coloured green) which are greater than two by one. Finally, we decrease any range labels larger than 4 (coloured blue) by one, to end up with a permutation in :math:`\mathcal{S}_5`.
		
		.. math::
			
			\begin{pmatrix}
				1	&\color{red}{2}	&3	&4	&5	&6	\\
				3	&\color{red}{4}	&5	&6	&1	&2
			\end{pmatrix}
			\to
			\begin{pmatrix}
				1	&3					&4					&5					&6	\\
				3	&\color{blue}{5}	&\color{blue}{6}	&1					&2
			\end{pmatrix}
			\\ \to 
			\begin{pmatrix}
				1	&\color{LimeGreen}{2}	&\color{LimeGreen}{3}	&\color{LimeGreen}{4}	&\color{LimeGreen}{5}	\\
				3	&\color{blue}{4}		&\color{blue}{5}		&1						&2
			\end{pmatrix}
		
		The same computation, but done in the interpreter:
		
			>>> x = Permutation("3 4 5 6 1 2"); print(x)
			(1 3 5)(2 4 6)
			>>> x.remove_image_of(2); print(x)
			(1 3 5 2 4)
			>>> x == Permutation("3 4 5 1 2")
			True
		
		Some further examples:
		
			>>> y = Permutation("2 5 3 1 4"); print(y)
			(1 2 5 4)
			>>> y.remove_image_of(2); print(y)
			(1 2 3)
			>>> z = Permutation("1 2 3 4 5 6") #identity of S_6
			>>> z.remove_image_of(4); print(z) #now identity of S_5
			()
		
		This method is mainly included for use by :meth:`~thompson.tree_pair.TreePair.reduce`.
		"""
		removed = self.output.pop(index-1)
		self.size -= 1
		for j, value in self:
			if value > removed:
				self[j] = value - 1
	
	def expand(self, index, to_width):
		r"""Expands the rule that index -> self[index] into *to_width* different rules: index -> self[index], index + 1 -> self[index] + 1, ...; then adjusts the labels used by the permutation to ensure it is still a bijection. This is used by :meth:`~thompson.tree_pair.TreePair.expand` when expanding a leaf into a caret for multiplication.
		
		**Example**. Suppose we begin with ``Permutation("3 4 5 2 1")`` and suppose we expands the image of 2 to a width of 3.
		We replace the symbol 2 by three copies (colored red) and assign images to them, counting up from the image of 2 (coloured blue). The twos are relabelled by integers, and the other indices are increased by 2 to compensate. To finish, we increase all original images which are greater than 5 by 2 (those in green) so that we have a bijection of {1, ..., 7}.
		
		.. math::
			
			\begin{pmatrix}
				1	&2	&3	&4	&5	\\
				3	&4	&5	&2	&1	
			\end{pmatrix}
			\to
			\begin{pmatrix}
				1	&\color{red}{2_A}	&\color{red}{2_B}	&\color{red}{2_C}	&3	&4	&5	\\
				3	&\color{blue}{4}	&\color{blue}{5}	&\color{blue}{6}		&5	&2	&1	
			\end{pmatrix}
			\\ \to 
			\begin{pmatrix}
				1	&\color{red}{2}		&\color{red}{3}		&\color{red}{4}		&5						&6	&7	\\
				3	&\color{blue}{4}	&\color{blue}{5}	&\color{blue}{6}	&\color{LimeGreen}{7}	&2	&1	
			\end{pmatrix}
		
		In the interpreter:
		
			>>> x = Permutation("3 4 5 2 1")
			>>> x.expand(2, to_width=3); x
			Permutation('3 4 5 6 7 2 1')
			
		Other examples:
		
			>>> y = Permutation("1 2 3 4");
			>>> y.expand(3, to_width=5); y #Expanding the identity gives the identity
			Permutation('1 2 3 4 5 6 7 8')
			>>> z = Permutation("5 3 2 6 1 4");
			>>> z.expand(3, to_width=3); z
			Permutation('7 5 2 3 4 8 1 6')
			>>> #expanding any index to a width of 1 does nothing
			>>> z.expand(3, to_width=1); z
			Permutation('7 5 2 3 4 8 1 6')
		"""
		if to_width < 1:
			raise ValueError("Width argument must be a positive integer (received {}).".format(to_width))
		image = self[index]
		
		#1. Increase the output values to make room in the range.
		for i, value in enumerate(self.output):
			if i + 1 != index and value > image:
				self.output[i] = value + to_width - 1
		
		replace = list(range(image, image + to_width))
		self.size += (to_width - 1)
		self.output = self.output[:index-1] + replace + self.output[index:] 
	

def lcm2(a, b):
	return a * b // gcd(a, b)

def lcm(iterable):
	"""Returns the least common multiple of *iterable*'s elements. From `Stack Overflow <http://stackoverflow.com/a/147539>`."""
	return reduce(lcm2, iterable, 1)
