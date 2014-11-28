from test import setup_script
setup_script(__file__)

import timeit

"""Repeatedly compute the QNBs of random automorphisms to see the average performance."""

num_examples = 10**4
repeats = 10


if __name__ == '__main__':
	runtime = timeit.timeit(
		'random_automorphism().quasinormal_basis()',
		'from thompson.examples import random_automorphism',
		number=num_examples)
	print('{} QNBs took {} seconds to compute. Average: {:.3f} ms.'.format(
		num_examples, runtime, 1000*runtime/num_examples))
	
	runtime = timeit.timeit(
		'phi.quasinormal_basis()',
		'from thompson import Automorphism; phi = Automorphism.from_file("phi_slow_QNB.aut")',
		number=1)
	print('Slow phi example took {:.2f} seconds to compute.'.format(runtime))