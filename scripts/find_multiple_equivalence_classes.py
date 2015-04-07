from test import setup_script
setup_script(__file__)

"""A script which seeks infinite automorphisms with more than one non-characteristic equivalence class."""

from random import randint

from networkx.algorithms.components import number_weakly_connected_components as num_components, \
                                           weakly_connected_components        as components

from thompson.examples import *

aut = load_example('multiple_classes_smaller')
print(aut)
aut.dump_QNB()
type_b, type_c = aut._split_basis()
roots, graph = aut.equivalence_data(type_b, type_c)
for i, comp in enumerate(components(graph)):
	print('component', i)
	for word in comp:
		print('\t', word)

def search():
	num_examples = 1
	while True:
		if num_examples % 100 == 0:
			print(num_examples)
		aut = random_infinite_automorphism(num_expansions=randint(1, 8))
		type_b, type_c = aut._split_basis()
		roots, graph = aut.equivalence_data(type_b, type_c)
		n = num_components(graph)
		if n == 1:
			num_examples += 1
			continue
		print("Found an example with {} equivalence classes".format(n))
		print(aut)
		aut.save_to_file('multiple_classes.aut')
		break
