from scripts import setup_script
setup_script(__file__)

from thompson.examples import *

"""Computes the products phi * phi^-1 and phi^-1 * phi and checks that they're both equal to the identity."""

i = 0
while True:
	i += 1
	if i % 100 == 0:
		print(i)
	phi = random_automorphism()
	e = ~phi * phi
	f = phi * ~phi
	if not (e == f and e.is_identity()):
		print(psi)
		break
