"""
.. testsetup::
	
	from thompson.periodic import *
	from thompson.examples import *
"""


from collections import defaultdict, deque

from .number_theory import lcm
from .generators    import Generators
from .automorphism  import Automorphism

__all__ = ["PeriodicAut"]

class PeriodicAut(Automorphism):
	r"""A purely periodic automorphism, which may have been extracted from a mixed automorphism.
	
		>>> print(example_5_9)
		PeriodicAut: V(2, 1) -> V(2, 1) specified by 7 generators (after expansion and reduction).
		x1 a1 a1 a1 -> x1 a1 a1 a2
		x1 a1 a1 a2 -> x1 a1 a2   
		x1 a1 a2    -> x1 a1 a1 a1
		x1 a2 a1 a1 -> x1 a2 a1 a2
		x1 a2 a1 a2 -> x1 a2 a1 a1
		x1 a2 a2 a1 -> x1 a2 a2 a2
		x1 a2 a2 a2 -> x1 a2 a2 a1
		>>> sorted(example_5_9.cycle_type)
		[2, 3]
		>>> #Two orbits of size 2, one orbit of size 3
		>>> from pprint import pprint
		>>> pprint(example_5_9.multiplicity)
		{2: 2, 3: 1}
		
	
	:ivar multiplicity: a mapping :math:`d \mapsto m_\phi(d, X_\phi)` where :math:`\phi` is the current automorphism and :math:`X_\phi` is the :meth:`quasi-normal basis <thompson.mixed.MixedAut.quasinormal_basis>` for :math:`\phi`.
	:ivar cycle_type: the set :math:`\{d \in \mathbb{N} : \text{$\exists$ an orbit of length $d$.}\}`
	:ivar order: the smallest positive number :math:`n` for which :math:`phi^n` is the identity. (This is the lcm of the cycle type.)
	"""
	def setup(self):
		assert self.quasinormal_basis.minimal_expansion_for(self) == self.quasinormal_basis
		#see definition 5.8
		self.cycle_type = set()
		counts = defaultdict(int)
		basis = self.quasinormal_basis
		
		for gen in basis:
			type, _, _ = self.orbit_type(gen, basis)
			orbit_size = type.characteristic[0]
			self.cycle_type.add(orbit_size)
			counts[orbit_size] += 1
		
		self.order = lcm(self.cycle_type)
		
		for orbit_size, occurances in counts.items():
			assert occurances % orbit_size == 0
			counts[orbit_size] = occurances // orbit_size
		self.multiplicity = dict(counts)
	
	def enumerate_orbits(self, basis):
		r"""Enumerates the periodic orbits of the current automorphism's quasinormal_basis. Returns a dictionary *orbits_by_size*. Each value ``orbits_by_size[d]`` is a list of the orbits of size *d*. Orbits themselves are represented as lists of :class:`Words <thompson.word.Word>`.
		
			>>> def display_orbits(orbits_by_size):
			... 	for key in sorted(orbits_by_size):
			... 		print('Orbits of length', key)
			... 		for list in orbits_by_size[key]:
			... 			print('...', *list, sep=' -> ', end=' -> ...\n')
			>>> orbits_by_size = example_5_9.enumerate_orbits(example_5_9.quasinormal_basis)
			>>> display_orbits(orbits_by_size)
			Orbits of length 2
			... -> x1 a2 a1 a1 -> x1 a2 a1 a2 -> ...
			... -> x1 a2 a2 a1 -> x1 a2 a2 a2 -> ...
			Orbits of length 3
			... -> x1 a1 a1 a1 -> x1 a1 a1 a2 -> x1 a1 a2 -> ...
		"""
		#TODO another doctest
		orbits_by_size = defaultdict(deque)
		already_seen = set()
		for gen in basis:
			if gen in already_seen:
				continue
			type, images, _ = self.orbit_type(gen, basis)
			length = type.characteristic[0]
			images = [images[i] for i in range(length)]
			already_seen.update(images)
			orbits_by_size[length].append(images)
		return dict(orbits_by_size)
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely periodic automorphisms are conjugate by examining their orbits.
		
			>>> psi_p = example_5_12_psi; phi_p = example_5_12_phi
			>>> rho_p = psi_p.test_conjugate_to(phi_p)
			>>> print(rho_p)
			PeriodicAut: V(2, 1) -> V(2, 1) specified by 6 generators (after expansion and reduction).
			x1 a1 a1 a1 a1 -> x1 a1 a2   
			x1 a1 a1 a1 a2 -> x1 a2 a2   
			x1 a1 a1 a2    -> x1 a1 a1 a1
			x1 a1 a2       -> x1 a2 a1 a1
			x1 a2 a1       -> x1 a1 a1 a2
			x1 a2 a2       -> x1 a2 a1 a2
			>>> rho_p * phi_p == psi_p * rho_p
			True
			
			>>> psi_p, phi_p = random_conjugate_periodic_pair()
			>>> rho_p = psi_p.test_conjugate_to(phi_p)
			>>> rho_p * phi_p == psi_p * rho_p
			True
		
		.. seealso:: This implements algorithm 5.13 of the paper---see section 5.3.
		"""
		if not isinstance(other, PeriodicAut):
			return None
		# 1. The quasi-normal bases are constructed in initialisation.
		# 2. Check that the cycle types are the same.
		if self.cycle_type != other.cycle_type:
			return None
		
		#3. Check that the multiplicities are congruent.
		modulus = self.signature.arity - 1
		for d in self.cycle_type:
			if self.multiplicity[d] % modulus != other.multiplicity[d] % modulus:
				return None
		
		# 4. Expand bases until the orbits multiplicites are the same
		s_orbits_of_size = self.enumerate_orbits(self.quasinormal_basis)
		o_orbits_of_size = other.enumerate_orbits(other.quasinormal_basis)
		
		for d in self.cycle_type:
			expansions_needed = (self.multiplicity[d] - other.multiplicity[d]) // modulus
			
			if expansions_needed > 0:
				expand_orbits(o_orbits_of_size[d], expansions_needed)
			elif expansions_needed < 0:
				expand_orbits(s_orbits_of_size[d], -expansions_needed)
			
			assert len(s_orbits_of_size[d]) == len(o_orbits_of_size[d])
		
		domain = Generators(self.signature)
		range  = Generators(self.signature)
		for d in self.cycle_type:
			for s_orbit, o_orbit in zip(s_orbits_of_size[d], o_orbits_of_size[d]):
				for s_word, o_word in zip(s_orbit, o_orbit):
					domain.append(s_word)
					range.append(o_word)
		
		rho = Automorphism(domain, range)
		rho.add_relabellers(self.domain_relabeller, other.range_relabeller)
		return rho
	
	def test_power_conjugate_to(self, other, multiple_solns=False):
		"""Tests two periodic factors to see if they are power conjugate. If *multiple_solns* is true, yields minimal soln (a, b, rho). Otherwise yields a single soln if it exists; otherwise None."""
		#This is almost exactly the same code as InfiniteAut.test_power_conjugate_to(). Maybe this should be one method on Automorphism
		
		if not isinstance(other, PeriodicAut):
			return None
		if self.is_identity() or other.is_identity():
			return None
		bounds = self.power_conjugacy_bounds(other)
		soln_iterator = self._test_power_conjugate_upto(other, *bounds, inverses=False)
		
		if multiple_solns:
			return soln_iterator
		else:
			try:
				return next(soln_iterator)
			except StopIteration:
				return None
	
	def power_conjugacy_bounds(self, other):
		"""We simply try all powers of both automorphisms. There are only finitely many, because everything is periodic.
		
		.. seealso:: Section 6.1."""
		return self.order, other.order

def expand_orbits(deque, num_expansions):
	r"""Takes a *deque* whose elements are a list of words. The following process is repeated *num_expansions* times.
	
	1. Pop an orbit :math:`\mathcal{O} = [o_1, \dotsc, o_d]` from the left of the deque.
	2. For :math:`i = 1, \dotsc, n` (where :math:`n` is the arity):
		
		a. Compute the new orbit :math:`\mathcal{O_i} =[o_1\alpha_i, \dotsc, o_d\alpha_i]`
		b. Append this orbit to the right of *deque*.
	
	Once complete, the number of orbits in *deque* has increased by :math:`\text{(arity -1) $\times$ num_expansions}`
	"""
	for _ in range(num_expansions):
		orbit = deque.popleft()
		arity = orbit[0].signature.arity
		for i in range(1, arity+1):
			new_orbit = [w.alpha(i) for w in orbit]
			deque.append(new_orbit)
