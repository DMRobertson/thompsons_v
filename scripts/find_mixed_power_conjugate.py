from test import setup_script
setup_script(__file__)

from thompson import Automorphism
# from thompson.mixed import MixedAut
# from thompson.examples import random_power_conjugate_pair

"""A script which seeks mixed power conjugate automorphisms"""

psi = Automorphism.from_file('mixed_pconj_psi.aut')
phi = Automorphism.from_file('mixed_pconj_phi.aut')

print(psi)
print(phi)
#orders are 3, 3
#power bounds for infinite are 625, 16
rho = psi.test_power_conjugate_to(phi)

def search():
	while True:
		psi, phi = random_power_conjugate_pair(signature=(2,1), num_expansions=5)
		print(psi.__class__.__name__, phi.__class__.__name__)
		if not (isinstance(psi, MixedAut) and isinstance(phi, MixedAut)):
			continue
		print(psi, phi)
		answer = input('continue?')
		if answer != 'y':
			break