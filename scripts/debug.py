from test import setup_script
setup_script(__file__)

from thompson.automorphism import Automorphism
from thompson.orbits import *
from thompson.examples import *
from pprint import pprint

from random import randint

psi = Automorphism.from_file('psi_tricky.aut')
print(repr(psi._qnf_basis))
print(psi)
basis = psi.quasinormal_basis()
print(basis)
for gen in basis:
	type, images, _ = psi.orbit_type(gen, basis)
	print(gen, type)
	# pprint(images)
p, i = psi._partition_basis(basis)
print(p)
print(i)