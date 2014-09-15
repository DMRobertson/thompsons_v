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
			output = [int(x) for x in output.split()]
		
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
			ValueError: "Key 10 is out of range 1, ..., 4."
		"""
		if not (1 <= key <= self.size):
			raise ValueError("Key %i is out of range 1, ..., %i." % (key, self.size))
		return self.output[key-1]
	
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
	
	def is_identity(self):
		"""Returns ``True`` if the permutation maps *x* to itself, for all applicable *x*. Otherwise ``False``.
		
			>>> Permutation("1 2 3 4").is_identity()
			True
			>>> Permutation("2 1 3 4").is_identity()
			False
		"""
		return all(self[i] == i for i in range(1, self.size+1))
	
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
	
	def remove_from_domain(self, index):
		r"""Removes the rule that index -> self[index] from the permutation. The permutation is relabelled so that it uses the symbols 1 to n.
		
		**Example**. Suppose we begin with ``Permutation("3 4 5 6 1 2")``. Removing 4 would give us the permutation described by the first matrix below. Next, we alter the map such that it sends the preimage of 4 to the image of 4 (coloured red). Then we decrease any symbols larger than 4 (coloured blue) by one, to end up with a permutation in :math:`\mathcal{S}_5`.
		
		.. math::
			
			\begin{pmatrix}
				1	&\color{red}{2}	&3	&\color{red}{-}	&5	&6	\\
				3	&\color{red}{-}	&5	&\color{red}{6}	&2	&1
			\end{pmatrix}
			\to
			\begin{pmatrix}
				1	&\color{red}{2}	&3	&-	&5	&6	\\
				3	&\color{red}{6}	&5	&-	&1	&2
			\end{pmatrix}
			\\\to 
			\begin{pmatrix}
				1	&2					&3					&\color{blue}{4}	&\color{blue}{5}	\\
				3	&\color{blue}{5}	&\color{blue}{4}	&1					&2
			\end{pmatrix}
			\hspace{1.7em}
		
		The same computation, but done in the interpreter:
		
			>>> x = Permutation("3 4 5 6 1 2"); print(x)
			(1 3 5)(2 4 6)
			>>> x.remove_from_domain(4); print(x)
			(1 3 4)(2 5)
			>>> y = Permutation("3 5 4 1 2"); print(y)
			(1 3 4)(2 5)
		
		This method is mainly included for use by :meth:`~thompson.tree_pair.TreePair.reduce`.
		"""
		#1. What maps to the thing we're going to remove?
		preimage = self.output.index(index) + 1
		#2. Remove the thing we want and fill in the hole we just created
		image = self.output.pop(index-1)
		self.size -= 1
		self.output[preimage-1] = image
		#3. Relabel so that everything fits from 1...n.
		for i, value in enumerate(self.output):
			if value > index:
				self.output[i] -= 1

def lcm2(a, b):
	return a * b // gcd(a, b)

def lcm(iterable):
	"""Returns the least common multiple of its arguments. From `Stack Overflow <http://stackoverflow.com/a/147539>`."""
	return reduce(lcm2, iterable, 1)
