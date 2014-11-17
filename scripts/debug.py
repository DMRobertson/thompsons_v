from test import setup_script
setup_script(__file__)

from thompson.automorphism import Automorphism
from thompson.orbits import *
from thompson.examples import *
from thompson.word import Word
from pprint import pprint

from random import randint

psi_i = Automorphism.from_file('psi_i.aut')
phi_i = Automorphism.from_file('phi_i.aut')

# print(psi_i)
# print(psi_i.quasinormal_basis())
# print(phi_i)
# print(phi_i.quasinormal_basis())

# rho_i = psi_i.test_conjugate_to(phi_i)

u = Word('x a1 a1 a1 a1 a1 a1 a1 a2', (2, 1))
print(phi_i.repeated_image(u, 4))
v = Word('x a2 a2 a1 a1 a2', (2, 1))
solns = phi_i.share_orbit(u, v)

print(solns)