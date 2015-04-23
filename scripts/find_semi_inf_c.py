from scripts import setup_script, find_examples_passing
setup_script(__file__)

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

def generate_automorphism():
	return random_automorphism(signature = (2, 1))

if __name__ == "__main__":
	find_examples_passing(test_semi_inf_c, generate_automorphism, 'semi_inf_c', description = "This test records examples of automorphisms whose quasinormal bases have elements of type C belonging to semi-infinite characteristics. (These are neccesary for ponds to exist but [I think] not neccesarily sufficient.)")