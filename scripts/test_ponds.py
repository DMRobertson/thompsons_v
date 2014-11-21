from test import setup_script
setup_script(__file__)

from thompson.word import Signature
from thompson.automorphism import Automorphism
from thompson.orbits import *
from thompson.examples import *
from pprint import pprint

import os
from random import randint


def dump_to_file(aut, name, log):
	print('{}\n====\n\n{}'.format(name.upper(), aut), file=log)
	if aut is None:
		try:
			os.remove(name + '.aut')
		except OSError:
			pass
	else:
		aut.save_to_file(name + '.aut')

num_attempts = 0
num_found = 0
exit = False
while True:
	num_attempts += 1
	if num_attempts % 100 == 0:
		print('tried', num_attempts, 'automorphisms')
	
	psi = random_automorphism()
	X = psi.quasinormal_basis()
	
	for endpt in psi.semi_infinite_end_points():
		type, _, _ = psi.orbit_type(endpt, X)
		if not type.is_type_C():
			continue
			
		LSI = type.type == OrbitType._left_semi_infinite
		sign = 1 if LSI else -1
		images = {i: psi.repeated_image(endpt, i) for i in range(-10, 11)}
		
		for i in range(2, 10):
			if X.is_above(images[i*sign]):
				print('found something')
				num_found += 1
				filename = 'pond_example_{}.txt'.format(num_found)
				
				with open(filename, 'wt') as f:
					print('details written to', filename)
					print(psi, file=f)
					pprint(images, f)
					print(endpt, type, file=f)
					print('Look at the power', i*sign, file=f)
				break
		if exit: break
	if exit: break
