from scripts import setup_script, find_examples_passing, choose_from_enum
setup_script(__file__)

import os
from enum import Enum

from thompson          import Automorphism
from thompson.periodic import PeriodicAut
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
	if aut.is_revealing(domain='wrt QNB'):
		return 'The QNB corresponds to a revealing pair.'
	return False

def is_not_revealing(aut):
	return not is_revealing(aut)

def nasty_snf_or_not_fully_informed(aut):
	X = aut.quasinormal_basis
	term, init = aut.semi_infinite_end_points()
	
	for x in term + init:
		if x not in X:
			continue
		ctype = aut.orbit_type(x)[0]
		if not ctype.is_type_B():
			#we have a SINC component in X
			return "Not nice: {} {}".format(x, ctype)
		power, mult = ctype.characteristic
		descendant = aut.repeated_image(x, power)
		assert x.is_above(descendant)
		
		domain = X.minimal_expansion_for(aut)
		range = aut.image_of_set(domain)
		
		leaves = range if power < 0 else domain
		backward = power < 0
		assert x in leaves, (x, ctype)
		iac_end = aut._end_of_iac(x, leaves, backward)
		if iac_end == descendant:
			continue
		#else:
		return "The component containing {} is not fully informed (descendant = {}, iac end = {}).".format(
		  x, descendant, iac_end)
	return False

def nice_snf_and_fully_informed(aut):
	return not nasty_snf_or_not_fully_informed(aut)

not_periodic = lambda f: not isinstance(f, PeriodicAut)

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

def find_snf_fully_informed():
	find_examples_passing(
	  test_functions         = [nice_snf_not_fully_informed],
	  automorphism_generator = generate_automorphism,
	  test_name              = 'snf_fully_informed',
	  description            = """We think that an SNF is associated to a revealing pair iff it is `fully informed'. This function finds examples of such SNFs.""",
	  max_examples           = 100,
	  save_examples          = True)

def find_fully_informed_non_revealing():
	find_examples_passing(
	  test_functions         = [nice_snf_and_fully_informed, is_not_revealing],
	  automorphism_generator = generate_automorphism,
	  test_name              = 'fully_informed_non_revealing',
	  description            = """We think that an SNF is associated to a revealing pair iff it is `fully informed'. This function searches for examples of such SNFs with are NOT revealing.""",
	  max_examples           = 100,
	  save_examples          = True)

def find_revealing_not_fully_informed():
	find_examples_passing(
	  test_functions         = [is_revealing, nasty_snf_or_not_fully_informed],
	  automorphism_generator = generate_automorphism,
	  test_name              = 'revealing_not_fully_informed',
	  description            = """We think that an SNF is associated to a revealing pair iff it is `fully informed'. This function searches for counter examples: revealing pairs which are not fully informed as SNFs.""",
	  max_examples           = 100,
	  save_examples          = True)

choices = {
	'find nice snf not revealing'      : find_nice_snf_not_revealing,
	'find nasty snf and revealing'     : find_nasty_snf_and_revealing,
	'find snf fully informed'          : find_snf_fully_informed,
	'find fully informed non revealing': find_fully_informed_non_revealing,
	'find revealing not fully informed': find_revealing_not_fully_informed
}

if __name__ == "__main__":
	choose_from_enum(choices)()
