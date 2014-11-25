from test import setup_script
setup_script(__file__)

from thompson.word import Signature
from thompson.orbits import *
from thompson.examples import *
from pprint import pprint

import os
from random import randint


def dump_to_file(aut, name, log):
	print('\n{}\n====\n\n{}'.format(name.upper(), aut), file=log)
	if aut is None:
		try:
			os.remove(name + '.aut')
		except OSError:
			pass
	else:
		aut.save_to_file(name + '.aut')

i = 0
while True:
	i += 1
	if i % 100 == 0:
		print(i)
	
	psi, phi, psi_p, phi_p, psi_i, phi_i = None, None, None, None, None, None
	try:
		psi, phi = random_conjugate_pair(signature=Signature(2,1))
		psi_p, psi_i = psi.free_factors()
		phi_p, phi_i = phi.free_factors()
		
		assert (psi_p is None) == (phi_p is None)
		if psi_p is not None:
			rho_p = psi_p.test_conjugate_to(phi_p)
			if rho_p is None or psi_p * rho_p != rho_p * phi_p:
				raise RuntimeError('Conjugacy failed for periodic.')
		
		assert (psi_i is None) == (phi_i is None)
		if psi_i is not None:
			rho_i = psi_i.test_conjugate_to(phi_i)
			if rho_i is None or psi_i * rho_i != rho_i * phi_i:
				print(repr(rho_i))
				raise RuntimeError('Problem with infinite')
	except Exception as e:
		with open('results.txt', 'wt', encoding='utf-8') as f:
			from traceback import print_exception
			import sys
			print_exception(type(e), e, None, file=f)
			print_exception(type(e), e, None, file=sys.stdout)
			dump_to_file(psi, 'psi', f)
			dump_to_file(psi_p, 'psi_p', f)
			dump_to_file(psi_i, 'psi_i', f)
			dump_to_file(phi, 'phi', f)
			dump_to_file(phi_p, 'phi_p', f)
			dump_to_file(phi_i, 'phi_i', f)
		break