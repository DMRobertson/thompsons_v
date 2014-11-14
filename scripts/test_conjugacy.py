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
	try:
		psi, phi = random_conjugate_pair()
		psi_p, psi_i = psi.free_factors()
		phi_p, phi_i = phi.free_factors()
		rho_p = psi_p.test_conjugate_to(phi_p)
		if rho_p is None or psi_p * rho_p != rho_p * phi_p:
			raise RuntimeError('Conjugacy failed for periodic.')
		# rho_i = psi_i.test_conjugate_to(phi_i)
		# if rho_i is None or psi_i * rho_i != rho_i * phi_i:
			# raise RuntimeError('Problem with infinite')
	except Exception as e:
		with open('results.txt', 'wt', encoding='utf-8') as f:
			from traceback import print_exception
			print_exception(type(e), e, None, file=f)
			import sys
			print_exception(type(e), e, None, file=sys.stdout)
			print
			print('PSI\n===\n\n', psi, file=f)
			print('PHI\n===\n\n', phi, file=f)
		psi.save_to_file('psi.aut')
		phi.save_to_file('phi.aut')
		break
