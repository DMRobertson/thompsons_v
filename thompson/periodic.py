"""
.. testsetup::
	
	from thompson.periodic import *
	from thompson.examples import *
"""

from collections import defaultdict, deque

from .generators    import Generators
from .automorphism  import Automorphism

__all__ = ["PeriodicAut"]

class PeriodicAut(Automorphism):
	r"""A purely periodic automorphism, which may have been extracted from a mixed automorphism.
	
		>>> example_5_9 = load_example('example_5_9')
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
		>>> example_5_9.order
		6
		
		>>> load_example('cyclic_order_six').order
		6
		>>> phi = random_periodic_automorphism()
		>>> 1 <= phi.order < float('inf')
		True
		>>> len(phi.cycle_type) != 0
		True
		>>> len(phi.characteristics)
		0
	"""
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely periodic automorphisms are conjugate by examining their orbits.
		
			>>> psi_p, phi_p = load_example('example_5_12_psi'), load_example('example_5_12_phi')
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
		
		.. seealso:: This implements algorithm :paperref:`alg:periodic` of the paper---see section :paperref:`subsectionP`.
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
		s_orbits_of_size = {
			d: deque(orbits) for d, orbits in self.periodic_orbits.items()
		}
		o_orbits_of_size = {
			d: deque(orbits) for d, orbits in other.periodic_orbits.items()
		}
		
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
	
	def find_power_conjugators(self, other, identity_permitted=False, cheat=False):
		#This is almost exactly the same code as InfiniteAut.test_power_conjugate_to(). Maybe this should be one method on Automorphism
		if not isinstance(other, PeriodicAut):
			raise StopIteration
		if not identity_permitted and (self.is_identity() or other.is_identity()):
			print("One of the automorphisms is the identity")
			raise StopIteration
		if cheat:
			from .examples.random import random_power_bounds
			bounds = random_power_bounds
		else:
			bounds = self.power_conjugacy_bounds(other, identity_permitted)
		yield from self._test_power_conjugate_upto(other, *bounds, inverses=False)
	
	def test_power_conjugate_to(self, other, cheat=False):
		r"""Tests two periodic factors to see if they are power conjugate. Yields minimal solutions :math:`(a, b, \rho)`. such that :math:`\rho^{-1}\psi^a\rho = \phi^b`.
		
		.. seealso:: Section :paperref:`torPower` of the paper.
		"""
		try:
			return next(self.find_power_conjugators(other, cheat=cheat))
		except StopIteration:
			return None
	
	def power_conjugacy_bounds(self, other, identity_permitted):
		"""We simply try all powers of both automorphisms. There are only finitely many, because everything is periodic.
		
		:returns: ``self.power, other.power`` if the identity is permitted as a solution; otherwise ``self.power - 1, other.power - 1``.
		
		.. todo:: Maybe this should be 0 to order if the identity is permitted and 1 to order otherwise?
		
		.. seealso:: Section :paperref:`torPower` of the paper."""
		if identity_permitted:
			return self.order, other.order
		else:
			return self.order - 1, other.order - 1

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
			new_orbit = tuple(w.alpha(i) for w in orbit)
			deque.append(new_orbit)
