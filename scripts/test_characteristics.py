from scripts import setup_script
setup_script(__file__)

from thompson.examples import *

import os

"""Tests a lemma that conjugate auts always have the same characteristics"""

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
	
	try:
		psi, phi = random_conjugate_pair()
		assert psi.characteristics == phi.characteristics
	except BaseException as e:
		with open('results.txt', 'wt', encoding='utf-8') as f:
			dump_to_file(psi, 'psi', f)
			dump_to_file(phi, 'phi', f)
		raise e