from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.word import format
from thompson.examples import *
from thompson.orbits import *
from pprint import pprint

import networkx as nx
# import matplotlib.pyplot as plt


def fmt_triple(edge_data):
	return "({}, {}, {})".format(
		format(edge_data['start_tail']), edge_data['power'], format(edge_data['end_tail'])
	)

ex = example_4_25_i
# ex = example_5_3_i
# ex = example_5_26_psi_i

print(ex)
print('QN basis:', ex.quasinormal_basis())

triples = ex.equivalence_classes()
for i, (graph, root, type_c_data) in enumerate(triples):
	print('component', i, 'with root', root)
	print('type B elements:', *graph.nodes_iter())
	print('type C elements:', *type_c_data)
	for node in graph:
		print('\tEdges out of', node)
		for source, target in graph.out_edges_iter(node):
			data = graph[source][target]
			print('\t\tto', target, 'with data\n\t\t\t', fmt_triple(data))

# G = nx.union_all(components)

# edge_labels = { (u, v) : fmt_triple(data) for (u, v, data) in G.edges_iter(data=True) }

# pos = nx.spring_layout(G)
# nx.draw(G,
	# node_size = 3000,
	# pos = pos,
	# with_labels = True
# )
# nx.draw_networkx_edge_labels(G, pos, edge_labels)
# plt.show()