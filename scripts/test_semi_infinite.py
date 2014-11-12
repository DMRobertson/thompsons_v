from test import setup_script
setup_script(__file__)

from thompson.word import Signature
from thompson.orbits import *
from thompson.examples import *
from pprint import pprint

from random import randint

i = 0
while True:
	i += 1
	if i % 100 == 0:
		print(i)
	num_expansions = randint(1, 20)
	signature = Signature(2,1)
	psi = random_automorphism(num_expansions, signature)
	basis = psi.quasinormal_basis()
	for gen in basis:
		type, _, _ = psi.orbit_type(gen, basis)
		if type.type == OrbitType._left_semi_infinite or OrbitType._right_semi_infinite and type.is_type_C():
			with open('output.txt', 'wt') as f:
				print(gen, 'is type C in QNF basis\n', basis, file=f)
				print(psi, file=f)
			break

