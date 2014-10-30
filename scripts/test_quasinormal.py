from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *
from pprint import pprint


ex = example_4_25_i
G = ex.equivalence_graph()
for node in G:
	print('Edges out of', node)
	for source, target in G.out_edges_iter(node):
		print('\tto', target, 'with data')
		print('\t\t', end='')
		pprint(G.get_edge_data(source, target))