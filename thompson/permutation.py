"""
.. testsetup:: 
	
	from thompson.permutation import *
"""

from fractions import gcd
from functools import reduce

__all__ = ["Permutation"]

class Permutation:
	r"""A permutation of a finite set :math:`\{1, \dotsc, N\}` of integers.
	
	:ivar size: the size :math:`N` of the set being permuted.
	"""
	def __init__(self, output):
		r"""This method accepts the *output* :math:`[f(1), f(2), \dotsc, f(N)]` of a permutation, speficied as
		
		- a list of integers, or
		- a space-separated string of base 10 integers (which is converted into the above)
		
		and creates a corresponding Permutation object.
		
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
		#TODO. Might be more natural to use __call__?
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
	
	def is_cycle(self):
		"""Returns ``True`` if the permutation consists of a single cycle, otherwise ``False``. 
		
			>>> Permutation("1 2 3 4").is_cycle() #identity
			True
			>>> Permutation("2 1 3 4").is_cycle() #2-cycle (1 2)
			True
			>>> Permutation("2 1 4 3").is_cycle() #product of 2-cycles (1 2)(3 4)
			False
		"""
		#Identity permutations' ((1)(2)...(N)) cycles() methods will return empty lists.
		return len(self.cycles()) < 2
	
	def cycles(self):
		r"""Returns a list of tuples representing this permutation as a product of cycles. "1-cycles" like (3) correspond to the identity permutation, and are not included in the output.
		
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
		"""Returns the order of the permutation. This is the smallest positive integer *n* such that the *n* th power of the permutation is the identity.
		
			>>> x = Permutation("2 3 1 4 6 5")
			>>> print(x)
			(1 2 3)(5 6)
			>>> x.order()
			6
		"""
		return lcm(len(cycle) for cycle in self.cycles())

def lcm2(a, b):
	return a * b // gcd(a, b)

def lcm(iterable):
	"""Returns the least common multiple of its arguments. From `Stack Overflow <http://stackoverflow.com/a/147539>`."""
	return reduce(lcm2, iterable, 1)