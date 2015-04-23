from scripts import setup_script, find_examples_passing
setup_script(__file__)

import os
from thompson          import Automorphism
from thompson.examples import random_automorphism

"""A script that hunts for ponds."""

def test_semi_inf_c(aut):
	X = aut.quasinormal_basis
	for x in X:
		ctype = aut.orbit_type(x, X)[0]
		if ctype.is_type_C() and ctype.type in {3, 4}:
			break
	else: #if we didn't break
		return ''
	return '\t {} belongs to a {}.'.format(x, ctype) 

def test_revealing(aut):
	domain = aut.quasinormal_basis.minimal_expansion_for(aut)
	if aut.is_revealing(domain):
		return 'The QNB DOES correspond to a revealing pair.'
	return ''

def generate_automorphism():
	return random_automorphism(signature = (2, 1))

def search():
	find_examples_passing(
	  test_functions         = (test_semi_inf_c, test_revealing),
	  automorphism_generator = generate_automorphism,
	  test_name              = 'semi_inf_c',
	  description            = """This test records examples of automorphisms whose quasinormal bases have elements of type C belonging to semi-infinite characteristics. (These are neccesary for ponds to exist but [I think] not neccesarily sufficient.)""",
	  max_examples           = 1)

if __name__ == "__main__":
	search()