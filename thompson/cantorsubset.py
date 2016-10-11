"""
.. testsetup::

	from thompson.cantorsubset    import *
	from thompson.examples.random import random_basis
"""

from .word       import Word, format_cantor, lambda_arguments
from .generators import Generators

__all__ = ["CantorSubset"]

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
		The order of the set is **not** neccesarily preserved, as the set is sorted by this function.
		
			>>> basis = random_basis(cls=CantorSubset)
			>>> basis.contract()
			>>> basis == CantorSubset.standard_basis(basis.signature)
			True
		"""
		self.sort()
		arity = self.signature.arity
		index = 0
		while index <= len(self) - arity:
			this = self[index]
			L = len(this)
			parent = Word(this[:L-1], self.signature, preprocess=False)
			if (
				L > 1
				and
				all( len(self[index+i]) == L for i in range(1, arity) )
				and
				all( self[index+i][:L-1] == parent for i in range(1, arity) )
			):
				self[index:index+arity] = [parent]
				#alpha_n is stored as -n
				index += 1 + min(0, parent[-1])
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
			[0, 1]
			>>> X = CantorSubset((2,1), ["1"])
			>>> X.simplify(); print(X)
			[1]
		"""
		
		self.expand_lambdas()
		self.contract()
		
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
	
	def __str__(self):
		r"""We use a notation introduced to us by Bleak: square brackets around a word stand for "the Cantor set underneath" its argument.
		We use the :func:`~thompson.word.format_cantor` function to display the elements of the generating set.
		
			>>> print(CantorSubset((2, 1), ["x a1"]))
			[0]
			>>> S = CantorSubset((2, 1), ["x a1", "x a1 a1", "x a2 a1 a1", "x a2 a1 a2"])
			>>> print(S)
			[0, 00, 100, 101]
			>>> S.simplify(); print(S)
			[0, 10]
		"""
		return "[" + ", ".join( format_cantor(word) for word in self ) + "]"
	
