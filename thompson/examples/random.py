from random import randint

from ..word import Signature
from ..generators import Generators
from ..automorphism import Automorphism

#todo docstrings
__all__ = ['random_signature', 'random_basis', 'random_automorphism', 'random_conjugate_pair']

def random_signature():
	return (randint(2, 4), randint(1, 5))

def random_basis(num_expansions, signature):
	basis = Generators.standard_basis(signature)
	
	for _ in range(num_expansions):
		i = randint(0, len(basis) - 1)
		basis.expand(i)
	return basis

def random_automorphism(num_expansions, signature):
	domain = random_basis(num_expansions, signature)
	range  = random_basis(num_expansions, signature)
	return Automorphism(domain, range)

def random_conjugate_pair(num_expansions=None, signature=None):
	if num_expansions is None:
		num_expansions = randint(1, 8)
	if signature is None:
		signature = random_signature()
	psi = random_automorphism(num_expansions, signature)
	rho = random_automorphism(num_expansions, signature)
	phi = ~rho * psi * rho
	return psi, phi