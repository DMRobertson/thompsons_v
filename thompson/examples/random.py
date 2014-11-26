from random import randint, shuffle 
from copy import copy

from ..word import Signature, Word
from ..generators import Generators
from ..automorphism import Automorphism
from ..factors import PeriodicFactor, InfiniteFactor

#todo docstrings
__all__ = ['random_signature', 'random_simple_word', 'random_basis',
	'random_automorphism', 'random_conjugate_pair',
	'random_conjugate_factors', 'random_conjugate_periodic_factors', 'random_conjugate_infinite_factors']

def needs_defaults(undec):
	def decd(num_expansions=None, signature=None):
		if num_expansions is None:
			num_expansions = randint(1, 5)
		if signature is None:
			signature = random_signature()
		return undec(num_expansions, signature)
	return decd

def random_signature():
	return Signature(randint(2, 4), randint(1, 5))

def random_simple_word(signature=None):
	signature = signature or random_signature()
	letters = []
	base = randint(1, signature.alphabet_size)
	letters.append(base)
	
	num_alphas = randint(0, 15)
	for i in range(num_alphas):
		letters.append(-randint(1, signature.arity))
	return Word(letters, signature)

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
	psi_p, psi_i = psi.free_factors()
	phi_p, phi_i = phi.free_factors()
	return psi_p, psi_i, phi_p, phi_i

@needs_defaults
def random_conjugate_periodic_factors(num_expansions, signature):
	psi_p = phi_p = None
	while psi_p is None or phi_p is None:
		psi_p, _, phi_p, _ = random_conjugate_factors()
	assert isinstance(psi_p, PeriodicFactor)
	return psi_p, phi_p

@needs_defaults
def random_conjugate_infinite_factors(num_expansions, signature):
	psi_i = phi_i = None
	while psi_i is None or phi_i is None:
		_, psi_i, _, phi_i = random_conjugate_factors()
	assert isinstance(psi_i, InfiniteFactor)
	return psi_i, phi_i