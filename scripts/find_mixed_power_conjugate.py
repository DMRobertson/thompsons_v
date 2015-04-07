from test import setup_script
setup_script(__file__)

from thompson.automorphism import Automorphism, PowerCollection
from thompson.mixed import MixedAut
from thompson.examples import *

"""A script which seeks mixed power conjugate automorphisms"""

psi, phi = load_example_pair('mixed_pconj')

print(psi)
print(phi)
assert psi.test_conjugate_to(phi) is None
assert (psi**-4).test_conjugate_to(phi**2) is not None
assert (psi**6).test_conjugate_to(phi**3) is not None
a, b, rho = psi.test_power_conjugate_to(phi)
print(a, b)
print(rho)

def search():
	while True:
		psi, phi = random_power_conjugate_pair(signature=(2,1), num_expansions=5)
		print(psi.__class__.__name__, phi.__class__.__name__)
		if not (isinstance(psi, MixedAut) and isinstance(phi, MixedAut)):
			continue
		print(psi, phi, sep='\n')
		answer = input('continue?')
		if answer != 'y':
			break