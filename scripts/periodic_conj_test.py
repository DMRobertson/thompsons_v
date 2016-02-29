from thompson          import Automorphism
from thompson.examples import random_conjugate_periodic_pair
import sys
from scripts import profile

def single_test():
	psi_p, phi_p = random_conjugate_periodic_pair()
	rho_p = psi_p.test_conjugate_to(phi_p)
	try:
		assert rho_p * phi_p == psi_p * rho_p
	except (KeyboardInterrupt, Exception) as e:
		psi_p.save_to_file('psi_p.aut')
		phi_p.save_to_file('phi_p.aut')
		raise e

def repeated_test():
	i = 0
	while True:
		i += 1
	#	if i % 50 == 0:
		print(i)
		single_test()

def debug():
	psi = Automorphism.from_file('psi_p.aut')
	phi = Automorphism.from_file('phi_p.aut')
	return psi, phi

if __name__ == '__main__':
	if 'debug' in sys.argv:
		psi, phi = debug()
	if 'profile' in sys.argv:
		psi, phi = debug()
		profile(psi.test_conjugate_to, args=[phi])
	else:
		repeated_test()
