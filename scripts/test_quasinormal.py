from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.word import format
from thompson.examples import *
from thompson.orbits import *
from pprint import pprint

def fmt_triple(edge_data):
	return "({}, {}, {})".format(
		format(edge_data['start_tail']), edge_data['power'], format(edge_data['end_tail'])
	)

def dump_graph(roots, graph):
	import networkx as nx
	for i, root in enumerate(roots):
		print('component', i, 'with root', root)
		print('type B elements:')
		for node in nx.dfs_preorder_nodes(graph, root):
			print('\tEdges out of', node)
			for source, target in graph.out_edges_iter(node):
				data = graph[source][target]
				print('\t\tto', target, 'with data\n\t\t\t', fmt_triple(data))

def plot_graph(G):
	import matplotlib.pyplot as plt
	edge_labels = { (u, v) : fmt_triple(data) for (u, v, data) in G.edges_iter(data=True) }
	pos = nx.spring_layout(G)
	nx.draw(G,
		node_size = 3000,
		pos = pos,
		with_labels = True
	)
	nx.draw_networkx_edge_labels(G, pos, edge_labels)
	plt.show()


psi = example_5_26_psi_i
phi = example_5_26_phi_i

# rho = psi.test_conjugate_to(phi)
print(phi.quasinormal_basis())
type_b, type_c = phi._split_basis()
roots, graph = phi.equivalence_data(type_b, type_c)
dump_graph(roots, graph)