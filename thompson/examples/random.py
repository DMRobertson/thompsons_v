from random import randint, shuffle
from copy import copy

from ..word import Signature
from ..generators import Generators
from ..automorphism import Automorphism
from ..factors import PeriodicFactor

#todo docstrings
__all__ = ['random_signature', 'random_basis',
	'random_automorphism', 'random_conjugate_pair',
	'random_conjugate_factors']

def needs_defaults(undec):
	def decd(num_expansions=None, signature=None):
		if num_expansions is None:
			num_expansions = randint(1, 5)
		if signature is None:
			signature = random_signature()
		return undec(num_expansions, signature)
	return decd

def random_signature():
	return (randint(2, 4), randint(1, 5))

@needs_defaults
def random_basis(num_expansions, signature):
	basis = Generators.standard_basis(signature)
	
	for _ in range(num_expansions):
		i = randint(0, len(basis) - 1)
		basis.expand(i)
	return basis

@needs_defaults
def random_automorphism(num_expansions, signature):
	domain = random_basis(num_expansions, signature)
	range  = random_basis(num_expansions, signature)
	shuffle(range)
	return Automorphism(domain, range)

@needs_defaults
def random_conjugate_pair(num_expansions, signature):
	psi = random_automorphism(num_expansions, signature)
	rho = random_automorphism(num_expansions, signature)
	phi = ~rho * psi * rho
	return psi, phi

@needs_defaults
def random_conjugate_factors(num_expansions, signature):
	psi, phi = random_conjugate_pair(num_expansions, signature)
	
