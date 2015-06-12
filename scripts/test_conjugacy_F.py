from scripts import setup_script, find_examples_passing
setup_script(__file__)

from functools         import lru_cache
from random            import randint
from thompson          import *
from thompson.examples import random_automorphism

"""Checks to see if our implementation of Higman's conjugacy test agrees with Hossain's implementation of Belk/Matucci's strand diagrams.
"""

signature = (2, 1)


def generate_auts():
	phi = random_automorphism(signature, group='F')
	psi = random_automorphism(signature, group='F')
	return phi, psi

# hossain_result = _____
# my_result = phi.test_conjugate_to(psi)
# return implies(hossain_result, my_result)

def implies(A, B):
	return B or not A
	
def test_via_hossain(phi, psi):
	phi = as_word(phi) 
	psi = as_word(psi) 
	#call his java applet here

def as_word(aut):
	"""See Cannon, Floyd, Parry 1996: Ex 2.4 and Thm 2.5"""
	domain_exponents = get_exponents(aut.domain)
	range_exponents  = get_exponents(aut.range)
	
	word = list(enumerate(range_exponents)) + list(enumerate(domain_exponents))[::-1]
	return ''.join(expand_letter(*letter) for letter in word)

@lru_cache(maxsize=None)
def expand_letter(sub, exp):
	letter = 'x' if exp > 0 else 'y'
	if sub in {0, 1}:
		return (letter + str(sub))*abs(exp)
	else:
		return expand_letter(0, 1-sub) + expand_letter(1, exp) + expand_letter(0, sub-1)

def get_exponents(basis):
	return [exponent(leaf, basis) for leaf in basis]

def exponent(leaf, basis):
	exponent = 0
	alphas = [ -x for x in leaf[1:] ]
	while len(alphas) > 1 and alphas.pop() == 1:
		if all(alpha == 2 for alpha in alphas):
			break
		exponent += 1
	return exponent

if __name__ == "__main__":
	phi = random_automorphism(signature=(2, 1), num_expansions=2, group="f")
	print(phi)
	print(as_word(phi))