from scripts import setup_script
setup_script(__file__)

import os
from pprint import pprint

from thompson.examples import *

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
	if True or i % 20 == 0:
		print('Test #', i, sep='')
	
	psi, phi = None, None
	try:
		print('Generating a pair')
		psi, phi = random_power_conjugate_pair()
		result = psi.test_power_conjugate_to(phi)
		assert result is not None
		pprint(result)
		a, b, rho = result
		assert (psi ** a) * rho == rho * (phi ** b)
	except BaseException as e:
		with open('results.txt', 'wt', encoding='utf-8') as f:
			f.write('Attempt number {}\n'.format(i))
			dump_to_file(psi, 'psi', f)
			dump_to_file(phi, 'phi', f)
		raise e