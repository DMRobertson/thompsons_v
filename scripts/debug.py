from test import setup_script
setup_script(__file__)

from thompson.automorphism import Automorphism
from thompson.orbits import *
from thompson.examples import *
from thompson.word import Word
from pprint import pprint
from itertools import chain

from random import randint

def dump_info():
	print('PSI')
	print(psi_i)
	print(psi_i.quasinormal_basis())
	print('PONDS: ', psi_i.pond_banks)
	print('PHI')
	print(phi_i)
	print(phi_i.quasinormal_basis())
	print('PONDS: ', phi_i.pond_banks)

psi_i = Automorphism.from_file('psi_i.aut')
phi_i = Automorphism.from_file('phi_i.aut')

dump_info()
dump_info()

print(psi_i.repeated_image(Word('x a2 a2 a1 a2 a2', (2,1)), 4))

rho_i = psi_i.test_conjugate_to(phi_i)
print(rho_i)
print(psi_i * rho_i == rho_i * phi_i)

dump_info()