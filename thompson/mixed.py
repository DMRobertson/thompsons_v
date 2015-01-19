"""
.. testsetup::
	
	from thompson.word import Word, Signature, from_string, free_monoid_on_alphas
	from thompson import word
	from thompson.generators import Generators
	from thompson.mixed import *
	from thompson.examples import *
	from thompson.number_theory import gcd
"""

__all__ = ["MixedAut"]

from copy import copy
from itertools import product, chain

from .word import *
from .generators import *
from .homomorphism import Homomorphism
from .automorphism import Automorphism
from .orbits import *

def modulo_non_zero(x, n):
	r"""Returns the unique integer :math:`s` such that :math:`1 \le s \le n` and :math:`s \equiv x \pmod{n}`.
	
		>>> [modulo_non_zero(x, 10) for x in range(20)]
		[10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
	"""
	x %= n
	if x == 0:
		return n
	return x

class MixedAut(Automorphism):
	#Extracting free factors
	def free_factors(self):
		p, i = self._partition_basis(self.quasinormal_basis)
		p_factor = self.free_factor(p, False)
		i_factor = self.free_factor(i, True)
		return p_factor, i_factor
	
	def _partition_basis(self, basis):
		r"""Let the current automorphism be in semi-normal form with respect to the given *basis*. This method places the elements of *basis* into two lists *(periodic, infinite)* depending on their orbit type. 
		
		:returns: the pair *(periodic, infinite)*. Both entries are sets of :class:`~thompson.generators.Generators`.
		
		.. doctest::
			
			>>> basis = example_4_5.quasinormal_basis
			>>> print(*example_4_5._partition_basis(basis), sep='\n')
			[x1 a2 a1, x1 a2 a2]
			[x1 a1 a1, x1 a1 a2]
			>>> basis = example_5_3.quasinormal_basis
			>>> print(*example_5_3._partition_basis(basis), sep='\n')
			[x1 a1 a2 a1, x1 a1 a2 a2]
			[x1 a1 a1 a1, x1 a1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> basis = example_5_3.quasinormal_basis
			>>> basis = basis.minimal_expansion_for(example_5_3)
			>>> print(*example_5_3._partition_basis(basis), sep='\n')
			[x1 a1 a2 a1, x1 a1 a2 a2]
			[x1 a1 a1 a1 a1, x1 a1 a1 a1 a2, x1 a1 a1 a2, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2]
		"""
		periodic = Generators(self.signature)
		infinite = Generators(self.signature)
		for gen in basis:
			type, _, _ = self.orbit_type(gen, basis)
			if type.is_incomplete():
				raise ValueError('MixedAut is not in semi-normal form with respect to the given basis.')
			elif type.is_type_A():
				periodic.append(gen)
			else:
				infinite.append(gen)
		return periodic, infinite
	
	def free_factor(self, generators, infinite=False):
		r"""This method restricts the current automorphism to the subalgebra generated by the given set :math:`X` of *generators*. This is then transformed into an automorphism of an isomorphic algebra with minimal alphabet size :math:`1 \le s \le n-1`.
		
		.. math::
		
			&G_{n, r}	&\to		&G_{n, |X|}		&\to		&G_{n, s} \\
			&\phi		&\mapsto	&\phi|_{\langle X\rangle}		&\mapsto	&\phi\,'
		
		:returns: The transformed automorphism :math:`\phi\, \in G_{n, s}`. Its type is :class:`PeriodicAut` if *infinite* is False; otherwise its type is :class:`InfiniteAut`.
		:raises ValueError: if an empty list of *generators* is provided.
		
		.. doctest::
			:hide:
			
			>>> qnb = example_5_3.quasinormal_basis
			>>> p, i = example_5_3._partition_basis(qnb)
			>>> print(example_5_3.free_factor(p, infinite=False))
			PeriodicAut: V(2, 1) -> V(2, 1) specified by 2 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 a2 a1 ~>    y1 a1 => y1 a2    ~> x1 a1 a2 a2
			x1 a1 a2 a2 ~>    y1 a2 => y1 a1    ~> x1 a1 a2 a1
		
		.. doctest::
			
			>>> # alphabet_size_two example
			>>> qnb = alphabet_size_two.quasinormal_basis
			>>> p, i = alphabet_size_two._partition_basis(qnb)
			>>> print(alphabet_size_two.free_factor(p, infinite=False))
			PeriodicAut: V(3, 1) -> V(3, 1) specified by 1 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 ~>    y1 => y1    ~> x1 a1
			>>> print(alphabet_size_two.free_factor(i, infinite=True))
			InfiniteAut: V(3, 1) -> V(3, 1) specified by 5 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a2 ~>    y1 a1    => y1 a1 a3    ~> x1 a2 a3
			x1 a3 ~>    y1 a2    => y1 a3       ~> x2      
			x2 a1 ~>    y1 a3 a1 => y1 a2       ~> x1 a3   
			x2 a2 ~>    y1 a3 a2 => y1 a1 a2    ~> x1 a2 a2
			x2 a3 ~>    y1 a3 a3 => y1 a1 a1    ~> x1 a2 a1
		"""
		if len(generators) == 0:
			raise ValueError('Must provide at least one generator.')
		#1. Decide how to relabel *generators* as elements of V_n,s
		modulus = self.signature.arity - 1
		alphabet_size = modulo_non_zero(len(generators), modulus)
		generators_relabelled = Generators.standard_basis((self.signature.arity, alphabet_size))
		generators_relabelled.expand_to_size(len(generators))
		
		#2. Relabel the domain and range
		domain = generators.minimal_expansion_for(self)
		range = self.image_of_set(domain)
		
		rewrite_set(domain, generators, generators_relabelled)
		rewrite_set(range,  generators, generators_relabelled)
		#Make sure we can undo the relabelling
		inverse_relabeller = Homomorphism(generators_relabelled, generators, reduce=False)
		
		#3. Return the factor
		factor = Automorphism(domain, range, reduce=False)
		factor.add_relabellers(inverse_relabeller, inverse_relabeller)
		#TODO We know a priori what the quasi-normal basis is. We could pass this info on to the factor?
		#TODO Would also be able to relabel type b data for bi-infinite orbit types
		return factor
	
	def _combine_factors(self, periodic, infinite):
		#doctest and docstring
		if periodic is None:
			p_domain = p_range = Generators(self.signature)
		else:
			p_domain, p_range = periodic.relabel()
		
		if infinite is None:
			i_domain = i_range = Generators(self.signature)
		else:
			i_domain, i_range = infinite.relabel()
		
		assert p_domain.signature == p_range.signature == i_domain.signature == i_range.signature == self.signature
		domain = p_domain + i_domain
		range  = p_range + i_range
		
		domain = Generators(self.signature, domain)
		range = Generators(self.signature, range)
		
		return MixedAut(domain, range)
	
	#Testing conjugacy
	def test_conjugate_to(self, other):
		r"""Determines if the current automorphism :math:`\psi` is conjugate to the *other* automorphism :math:`\phi`.
		
		:returns: if it exists, a conjugating element :math:`\rho` such that :math:`\rho^{-1}\psi\rho = \phi`. If no such :math:`\rho` exists, returns ``None``.
		:raises ValueError: if the automorphisms have different arities or alphabet sizes.
		
		.. doctest::
			
			>>> psi, phi = random_conjugate_pair()
			>>> rho = psi.test_conjugate_to(phi)
			>>> rho is not None
			True
			>>> psi * rho == rho * phi
			True
			>>> psi, phi = first_pond_example_psi, first_pond_example_phi
			>>> rho = psi.test_conjugate_to(phi)
			>>> rho is not None
			True
			>>> psi * rho == rho * phi
			True
		
		..doctest::
			:hide:
			
			>>> rho = nathan_pond_example.test_conjugate_to(nathan1_example)
			>>> rho is None
			True
			>>> rho = nathan2_example.test_conjugate_to(nathan3_example)
			>>> rho is not None
			True
			>>> nathan2_example * rho == rho * nathan3_example
			True
			>>> rho = nathan4_example
			>>> nathan2_example * rho == rho * nathan3_example
			True
			
		
		.. seealso:: This is an implementation of Algorithm 5.6 in the paper. It depends on Algorithms 5.13 and 5.27; these are the :meth:`periodic <thompson.factors.PeriodicAut.test_conjugate_to>` and :meth:`infinite <thompson.factors.InfiniteAut.test_conjugate_to>` tests for conjugacy.
		"""
		#TODO Doctest: try assembling a conjugator from factors
		#0. Check that both automorphisms belong to the same G_{n,r}.
		if self.signature != other.signature:
			raise ValueError('MixedAut\'s signatures {} and {} do not match.'.format(
			  self.signature, other.signature))
		
		#1. Before going any further, check that the number of periodic and infinite elements are compatible.
		sizes_okay = self._check_parition_sizes(other)
		if not sizes_okay:
			return None
		
		s_p, s_i = self.free_factors()
		o_p, o_i = other.free_factors()
		
		rho_p = s_p.test_conjugate_to(o_p)
		if rho_p is None:
			return None
		
		rho_i = s_i.test_conjugate_to(o_i)
		if rho_i is None:
			return None
		
		#Step 6. If we've got this far, we have a periodic/infinite mixture.
		#Combine the rewritten conjugators into one conjugator in G_n,r
		return self._combine_factors(rho_p, rho_i)
	
	def _check_parition_sizes(self, other):
		r"""Checks the sizes of a partition of two quasinormal bases to see if they might possibly describe conjugate automorphisms.
		
		:returns: None if we discover conjugacy is impossible. Otherwise, returns a pair of booleans *(pure_periodic, pure_infinite)*. Note that even if such a pair is returned, we do **not** know for certain that conjugacy is possible.
		"""
		s_qnf_p, s_qnf_i = self._partition_basis(self.quasinormal_basis)
		o_qnf_p, o_qnf_i = other._partition_basis(other.quasinormal_basis)
		
		#Check that the lengths match up modulo n-1.
		modulus = self.signature.arity - 1
		size_check = (len(s_qnf_p) % modulus == len(o_qnf_p) % modulus
		  and   len(s_qnf_i) % modulus == len(o_qnf_i) % modulus)
		
		return size_check
	
	#Power conjugacy
	def test_power_conjugate_to(self, other):
		r"""Determines if some power of the current automorphism :math:`\psi` is conjugate to some power of the *other* automorphism :math:`\phi`.
		
		:returns: if it exists, a triple :math:`(a, b, \rho)` such that :math:`\rho^{-1}\psi^a\rho = \phi^b`. If no such triple exists, returns ``None``.
		:raises ValueError: if the automorphisms have different arities or alphabet sizes.
		
		.. seealso:: This is an implementation of Algorithm 6.12 in the paper. It depends on Algorithms 5.6 (the :meth:`conjugacy test <test_conjugate_to>`) and 6.10 (the :meth:`infinite power conjugate test <thompson.factors.InfiniteAut.find_power_conjugators>.`)
		"""
		#0. Check that both automorphisms belong to the same G_{n,r}.
		if self.signature != other.signature:
			raise ValueError('MixedAut\'s signatures {} and {} do not match.'.format(
			  self.signature, other.signature))
		
		#1. Before going any further, check that the number of periodic and infinite elements are compatible.
		sizes_okay = self._check_parition_sizes(other)
		if not sizes_okay:
			return None
		
		#2. Periodic minimal solns.
		periodic_conjugators = list(s_p.test_power_conjugate_to(o_p, multiple_solns=True))
		if len(periodic_conjugators) == 0:
			return None
		
		#3. Infinite minimal solns.
		infinite_conjugators = list(s_i.test_power_conjugate_to(s_i, multiple_solns=True))
		if len(infinite_conjugators) == 0:
			return None
		
		#5. Try to recombine.
		for alpha, beta, rho_i in infinite_conjugators:
			for c, d, rho_p in periodic_conjugators:
				solns = solve_linear_congruence(a, c, s_p.order) & solve_linear_congruence(b, d, o_p.order)
				if not solns.is_empty():
					soln = solns.base
					rho = self._combine_factors(rho_p, rho_i)
					return alpha*soln, beta*soln, rho
		
		#6. If we've got this far, we're out of luck.
		return None
