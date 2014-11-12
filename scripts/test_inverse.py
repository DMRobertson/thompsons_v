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
	phi, _= random_conjugate_pair()
	e = ~phi * phi
	f = phi * ~phi
	if not (e == f and e.is_identity()):
		print(psi)
		break
