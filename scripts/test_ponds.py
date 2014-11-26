from test import setup_script
setup_script(__file__)

from thompson.orbits import *
from thompson.examples import *

from pprint import pprint
import os

"""A script that hunts for ponds"""

num_attempts = 0
num_found = 0
break_outer = False
while True:
	num_attempts += 1
	if num_attempts % 100 == 0:
		print('tried', num_attempts, 'automorphisms')	
	
	psi = random_automorphism()
	X = psi.quasinormal_basis()
	
	terminal, initial = psi.semi_infinite_end_points()
	for endpt in terminal:
		type, _, _ = psi.orbit_type(endpt, X)
		if not type.is_type_C():
			continue
			
		LSI = type.type == ComponentType._left_semi_infinite
		sign = 1 if LSI else -1
		images = {i: psi.repeated_image(endpt, i) for i in range(-10, 11)}
		
		for i in range(2, 10):
			if X.is_above(images[i*sign]):
				print('found something')
				num_found += 1
				filename = 'ponds/pond_example_{}.txt'.format(num_found)
				os.makedirs('ponds', exist_ok=True)
				
				with open(filename, 'wt') as f:
					print('details written to', filename)
					print(psi, file=f)
					pprint(images, f)
					print(endpt, type, file=f)
					print('Look at the power', i*sign, file=f)
				break_outer = True
				break
			if break_outer: break
		if break_outer: break
	break_outer = False

"""
The following regex will search for ponds of width (at least) 2

	L', \(\d, \d\)\),.*L', \(\d, \d\)\),

Note that '. matches newline' should be enabled in N++'s find dialogue.
"""