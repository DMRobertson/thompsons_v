"""TODO.

.. testsetup:: 
	
	from thompson.permutation import *
"""

class Permutation:
	r"""A permutation of the finite set :math:`\{1, \dotsc, N}` of integers.
	
	:ivar size: the size :math:`N` of the set being permuted.
	"""
	def __init__(self, output):
		"""Given the *output* of a permutation, described as either
		
		- a list of integers, or
		- a space-separated string of base 10 integers (which is converted into the above),
		
		creates a Permutation object mapping 1 to the first element of output[i], 2 to the second element, etc. Because Python uses zero-based indexing, ``i`` maps to ``output[i-1]``.
		
		:raises ValueError: if the given output do not properly define a permutation. For example:
		
			>>> Permutation('1 2 1 3')
			Traceback (most recent call last):
			...
			ValueError: Label 1 reused in ['1', '2', '1', '3']
			>>> Permutation('1 5 2 7')
			Traceback (most recent call last):
			...
			ValueError: Elements [3, 4] are not included in the output [1, 5, 2, 7]
			>>> Permutation('1 2 1 3 3')
			Traceback (most recent call last):
			...
			ValueError: Label '1' reused in [1 2 1 3 3]
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
	
	def __repr__(self):
		out_str = " ".join(str(x) for x in self.output)
		return self.__class__.__name__ + "(" + out_str + ")"
	
	def __getitem__(self, key):
		"""Permutations can be accessed like a dictionary to map inputs to output.
		
			>>> x = Permutation("4 3 1 2")
			>>> x[1]
			4
			>>> x[10]
			Traceback (most recent call last):
			...
			ValueError: "Key 10 is out of range 1, ..., 4."
		"""
		if not (1 <= key <= self.size):
			raise ValueError("Key %i is out of range 1, ..., %i." % (key, self.size))
		return self.output[key-1]

