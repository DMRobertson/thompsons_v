from test import setup_script
setup_script(__file__)

from pprint import pprint
import os

from thompson.examples import *

"""A script that hunts for ponds."""

num_attempts = 0
num_found = 0
break_outer = False
while True:
	num_attempts += 1
	if num_attempts % 10000 == 0:
		print('tried', num_attempts, 'automorphisms')
	
	psi = random_automorphism()
	X = psi.quasinormal_basis()
	
	if len(psi.pond_banks) == 0:
		continue
	
	print('Found {} ponds of lengths {} in an element of G_{}'.format(
		len(psi.pond_banks), {bank[1] - 1 for bank in psi.pond_banks}, psi.domain.signature))
	num_found += 1
	os.makedirs('ponds', exist_ok=True)
	filename = 'ponds/pond_example_{}.txt'.format(num_found)
	
	with open(filename, 'wt') as f:
		print('details written to', filename)
		print(psi, file=f)
		for (ell, k, r) in psi.pond_banks:
			print(*(psi.repeated_image(ell, i) for i in range(-1, k+2)), sep=', ', file=f)
"""
The following regex will search for ponds of width (at least) 2

	L', \(\d, \d\)\),.*L', \(\d, \d\)\),

Note that '. matches newline' should be enabled in N++'s find dialogue.
"""