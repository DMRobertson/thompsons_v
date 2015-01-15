from test import setup_script
setup_script(__file__)

from thompson.examples import random_automorphism
import time

"""Look for examples of QNBs that are slow to compute."""

num_examples = 10**6
num_repeats = 3
max_runtime = -1

slowest_aut = None

for i in range(1, num_examples+1):
	aut = random_automorphism()
	start = time.time()
	aut.quasinormal_basis
	runtime = time.time() - start
	
	if runtime > max_runtime:
		print('Slowest runtime thus far is {:2f}s'.format(runtime))
		max_runtime = runtime
		slowest_aut = aut
		aut.save_to_file('slow_example.aut')
	
	if i % 1000 == 0:
		print('tried {} examples'.format(i))

print(slowest_aut)