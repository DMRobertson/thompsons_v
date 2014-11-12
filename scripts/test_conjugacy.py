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
	psi, phi = random_conjugate_pair()
	rho = psi.test_conjugate_to(phi)
	if rho is None or psi * rho != rho * phi:
		print(psi, '\n\n\n', phi)
		break
