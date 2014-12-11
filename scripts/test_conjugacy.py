from test import setup_script
setup_script(__file__)

from thompson.examples import *

import os

"""Runs the conjugacy test on randomly generated conjugate pairs to ensure that the conjugacy test isn't producing false negatives."""

def dump_to_file(aut, name, log):
	print('\n{}\n====\n\n{}'.format(name.upper(), aut), file=log)
	if aut is None: #can pass in None if a periodic factor doesn't exist.
		try:
			os.remove(name + '.aut')
		except OSError:
			pass
	else:
		aut.save_to_file(name + '.aut')

i = 0
while True:
	i += 1
	if i % 20 == 0:
		print(i)
	
	psi, phi = None, None
	try:
		psi, phi = random_conjugate_pair()
		rho = psi.test_conjugate_to(phi)
		assert rho is not None
		assert psi * rho == rho * phi
	except BaseException as e:
		with open('results.txt', 'wt', encoding='utf-8') as f:
			dump_to_file(psi, 'psi', f)
			dump_to_file(phi, 'phi', f)
		raise e