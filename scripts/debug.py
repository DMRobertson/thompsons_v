from test import setup_script
setup_script(__file__)

from thompson.automorphism import Automorphism
from thompson.orbits import *
from thompson.examples import *
from thompson.word import Word
from pprint import pprint
from itertools import chain

from random import randint

psi_i = Automorphism.from_file('first_pond_example_psi.aut')
phi_i = Automorphism.from_file('first_pond_example_phi.aut')

print(psi_i)
print(psi_i.quasinormal_basis())
print(phi_i)
print(phi_i.quasinormal_basis())

# basis = phi_i.quasinormal_basis()
# for gen in chain(*phi_i.semi_infinite_end_points()):
	# type = phi_i.orbit_type(gen, basis)[0]
	# if not type.is_type_C():
		# continue
	# print('\n\n', gen, type)
	# for i in range(-10, 11):
		# print('\tpower', i, '=', phi_i.repeated_image(gen, i))

rho_i = psi_i.test_conjugate_to(phi_i)
print(rho_i)

# u = Word('x a1 a1 a1 a1 a1 a1 a1 a2', (2, 1))
# print(phi_i.repeated_image(u, 4))
# v = Word('x a2 a2 a1 a1 a2', (2, 1))
# solns = phi_i.share_orbit(u, v)

# print(solns)