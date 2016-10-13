"""
.. testsetup::

	from thompson.cantorsubset    import *
	from thompson.examples.random import random_basis, random_generators
	from thompson.word            import Word
	import random
"""

from .word       import Word, format_cantor, lambda_arguments
from .generators import Generators

__all__ = ["CantorSubset", "detailed_comparison"]

class CantorSubset(Generators):
	r"""A subclass of :class:`~thompson.generators.Generators` specifically designed to represent an open subset of the Cantor set :math:`\mathcal C`.
	We store this as a (typically sorted) list of :class:`~thompson.generators.Generators` :math:`u_1, \dots, u_n`.
	This stands for the set of all points :math:`u_i\Gamma` where :math:`\Gamma \in \{0,1\}^{\mathbb N}`.
	"""
	def expand_lambdas(self):
		"""
			>>> X = CantorSubset((2,1), ["x1 a1 a1 x1 a1 a2 a1 x1 a2 a1 a1 a1 L L"])
			>>> X.expand_lambdas(); print(X);
			[00, 010, 1000]
		"""
		index = 0
		while index < len(self):
			if not self[index].is_simple():
				self[index:index+1] = lambda_arguments(self[index])
			index += 1
	
	def contract(self):
		r"""Contracts the current generating set as much as possible (without using words involving a :math:`\lambda`).
		The set should be sorted before using this method.
		
			>>> basis = random_basis(cls=CantorSubset)
			>>> basis.contract()
			>>> basis == CantorSubset.standard_basis(basis.signature)
			True
		"""
		arity = self.signature.arity
		index = 0
		while index <= len(self) - arity:
			this = self[index]
			L = len(this)
			parent = Word(this[:L-1], self.signature, preprocess=False)
			if (
				all( len(self[index+i]) == L > 1 for i in range(1, arity) )
				and
				all( self[index+i][:L-1] == parent for i in range(1, arity) )
			):
				self[index:index+arity] = [parent]
				#alpha_n is stored as -n
				index += 1
				if len(parent) > 1:
					index += min(0, parent[-1])
			else:
				index += 1
	
	def simplify(self):
		r"""Simplifies the current Cantor subset to a normal form.
		We do this in three steps.
		Firstly we expand any :meth:`nonsimple Words <thompson.word.Word.is_simple>`.
		Second, we :meth:`contract` as much as possible (without using words involving a :math:`\lambda`).
		Last, we check to see if the set contains any pairs :math:`u, u\Gamma` and remove the latter.
		
			>>> X = CantorSubset((2,1), "0 11 10 00 111 1101 11110".split())
			>>> X.simplify(); print(X)
			[<entire Cantor set>]
			>>> X = CantorSubset((2,1), ["1"])
			>>> X.simplify(); print(X)
			[1]
		"""
		self.expand_lambdas()
		self.sort()
		
		index = 0
		while index < len(self) - 1:
			word = self[index]
			next_index = index + 1
			while next_index < len(self):
				next = self[next_index]
				if next[:len(word)] == word:
					next_index += 1
				else:
					break
			del self[ index+1:next_index ]
			index += 1
		
		self.contract()
	
	def __str__(self):
		r"""We use a notation introduced to us by Bleak: square brackets around a word stand for "the Cantor set underneath" its argument.
		We use the :func:`~thompson.word.format_cantor` function to display the elements of the generating set.
		
			>>> print(CantorSubset((2, 1)))
			[]
			>>> print(CantorSubset((2, 1), ["x"]))
			[<entire Cantor set>]
			>>> print(CantorSubset((2, 1), ["x a1"]))
			[0]
			>>> S = CantorSubset((2, 1), ["x a1", "x a1 a1", "x a2 a1 a1", "x a2 a1 a2"])
			>>> print(S)
			[0, 00, 100, 101]
			>>> S.simplify(); print(S)
			[0, 10]
		"""
		return "[" + ", ".join( format_cantor(word, subset=True) for word in self ) + "]"
	
	def __and__(self, other):
		r"""Computes the intersection of Cantor subsets.
		We assume the list of words is sorted before calling this function.
		
			>>> A = CantorSubset((2, 1), ["00"])
			>>> B = CantorSubset((2, 1), "01 111".split())
			>>> C = A & B
			>>> print(A, B, C)
			[00] [01, 111] []
			>>> D = CantorSubset((2, 1), ["11"])
			>>> print(D, B, D & B)
			[11] [01, 111] [111]
		
		"""
		if not isinstance(other, type(self)):
			raise TypeError("Don't know how to intersect {} and {} instances".format(
			  type(self).__name__, type(other).__name__))
		
		return type(self)(self.signature, ( x for x in self._intersect(other) ))
	
	def _intersect(self, other):
		sindex = 0
		oindex = 0
		while sindex < len(self) and oindex < len(other):
			subword, comparison = detailed_comparison(self[sindex], other[oindex]) 
			if subword:
				if comparison == 1:
					yield self[sindex]
					sindex += 1
					self, other    = other, self
					sindex, oindex = oindex, sindex
				elif comparison == -1:
					yield other[oindex]
					oindex += 1
				else:
					yield self[sindex]
					sindex += 1
					oindex += 1
			else:
				if comparison == 1:
					oindex += 1
				else:
					sindex += 1
	
	def __invert__(self):
		"""The complement. again it only works on a sorted list of leaves.
		
			>>> X = CantorSubset((2, 1), "01 1010 11".split())
			>>> print(~X)
			[00, 100, 1011]
			>>> print(~CantorSubset((2, 1), []))
			[<entire Cantor set>]
			>>> print(~CantorSubset((2, 1), ["x1"]))
			[]
			
			>>> X = random_generators(cls=CantorSubset, signature=(2,1))
			>>> X.sort();
			>>> comp = ~X
			>>> print(X & comp)
			[]
			>>> X.extend(comp); X.simplify()
			>>> print(X)
			[<entire Cantor set>]
		"""
		assert self.signature.alphabet_size == 1
		return type(self)(self.signature, ( Word(x, self.signature) for x in self._complement() ))
	 
	def _complement(self):
		def print(*args, **kwargs): pass
		if len(self) == 0:
			yield "x1"
			return
		
		needle = [1] #the root x1
		index = 0
		while index < len(self):
			target = self[index]
			subword, comparison = detailed_comparison(needle, target)
			print("needle:", format_cantor(needle), "target:", format_cantor(target), "subword:", subword, "comparison:", comparison)
			if subword and comparison == -1:
				#needle strictly above word
				print("attempting to add on zeroes")
				extend_zeroes(needle, target)
				print("extended needle:", format_cantor(needle), "target:", format_cantor(target))
				if len(target) == len(needle):
					print("extended needle same length as target")
					comparison = 0
				else:
					assert len(target) > len(needle)
					print("extended needle shorter than target, adding a zero")
					needle.append(-1)
					print("extra zero on needle:", format_cantor(needle))
					subword = False
					comparison = -1
				
			if subword and comparison == 0:
				needle == target
				print("needle == target\n===DON'T YIELD")
				needle = next_leaf(needle)
				index += 1
				print("next target, next needle")
			elif not subword and comparison == -1:
				print("needle not above target")
				print("===YIELD ", format_cantor(needle))
				yield needle
				needle = next_leaf(needle)
				print("next needle", format_cantor(needle))
			else:
				raise Exception("This shouldn't happen. Design your algorithm more carefully next time David!")
			print("bottom of loop")
		
		while needle is not None:
			yield needle
			needle = next_leaf(needle)

def extend_zeroes(needle, target):
	start = len(needle)
	index = start
	while index < len(target) and (target[index] == -1): #alpha 1
		needle.append(-1)
		index += 1
	return index - start

def next_leaf(letters):
	while True:
		if len(letters) == 1 and letters[0] > 0: #root
			return None
		elif letters[-1] == -1: # should be > -arity here
			letters[-1] = -2
			return letters
		else:
			del letters[-1]

def detailed_comparison(self, other):
	"""Returns a tuple *(subword, comparison)* of two Booleans.
	Either one word is a subword of the other or not.
	In the former case *subword* is ``True``, otherwise *subword* is ``False``.
	
	In both cases we can decide if *self* is lexicographically smaller than *other*.
	If so, *comparison* is ``-1``. If they are equal, *comparison* is ``0``; if *self* is larger than *other* then *comparison* is ``1``.
	
		>>> for thing in ["0 00", "001 00", "101 101", "0 1", "11 00"]:
		... 	thing = thing.split()
		... 	s = Word(thing[0], (2, 1))
		... 	o = Word(thing[1], (2, 1))
		... 	print(*detailed_comparison(s, o))
		... 
		True -1
		True 1
		True 0
		False -1
		False 1
	"""
	for s, o in zip(self, other):
		result = cmp(s, o)
		if result != 0:
			return ( False, -result )
	#If we reach here, the one is a subword of the other
	return ( True, cmp(len(self), len(other)) ) 

def cmp(a, b):
	return (a > b) - (a < b)
	