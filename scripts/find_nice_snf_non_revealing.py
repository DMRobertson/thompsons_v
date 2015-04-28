from scripts import setup_script, find_examples_passing, choose_from_enum
setup_script(__file__)

import os
from enum import Enum

from thompson          import Automorphism
from thompson.examples import random_automorphism

"""Not all semi-normal forms correspond to revealing pairs. This file contains scripts to try and explain which SNFs aren't revealing."""

def contains_semi_inf_c(aut):
	"""This test looks for auts for which no element x of the quasinormal_basis X belongs to a semi_infinite non-characteristic X-component."""
	X = aut.quasinormal_basis
	for x in X:
		ctype = aut.orbit_type(x, X)[0]
		if ctype.is_type_C() and ctype.type in {3, 4}:
			break
	else: #if we didn't break
		return False
	return '\t {} belongs to a {}.'.format(x, ctype) 

def no_semi_inf_c(aut):
	return not contains_semi_inf_c(aut)

def is_revealing(aut):
	"""This test looks for auts which are NOT revealing wrt their quasinormal_basis."""
	if aut.is_revealing(aut.quasinormal_basis.minimal_expansion_for(aut)):
		return 'The QNB corresponds to a revealing pair.'
	return False

def is_not_revealing(aut):
	return not is_revealing(aut)

def generate_automorphism():
	return random_automorphism(signature = (2, 1))

def find_nice_snf_not_revealing():
	find_examples_passing(
	  test_functions         = [no_semi_inf_c, is_not_revealing],
	  automorphism_generator = generate_automorphism,
	  test_name              = 'nice_snf_non_revealing',
	  description            = """We think that the SNF associated to a revealing pair contains no SI non-characteristic (sinc) elements. Does every SNF of this form yield a revealing pair? We test this by looking for counter-examples.""",
	  max_examples           = 100,
	  save_examples          = True)

def find_nasty_snf_and_revealing():
	find_examples_passing(
	  test_functions         = [contains_semi_inf_c, is_revealing],
	  automorphism_generator = generate_automorphism,
	  test_name              = 'nasty_snf_is_revealing',
	  description            = """We think that the SNF associated to a revealing pair contains no SI non-characteristic (sinc) elements. If the SNF is not of this form, can it ever be revealing? We search for examples.""",
	  max_examples           = 100,
	  save_examples          = True)

def examine(index):
	f = Automorphism.from_file('revealing/revealing_{}.aut'.format(index))
	print(f)
	f.dump_QNB()
	print('Is minimal rep revealing?', f.is_revealing())
	print('Evidence for non-revealing:', f.test_revealing())
	Y = f.quasinormal_basis.minimal_expansion_for(f)
	print('Is QNB revealing?', f.is_revealing(Y))
	print('Evidence for non-revealing:', f.test_revealing(Y))

choices = {
	'find nice snf not revealing' : find_nice_snf_not_revealing,
	'find nasty snf and revealing': find_nasty_snf_and_revealing,
}

if __name__ == "__main__":
	choose_from_enum(choices)()
