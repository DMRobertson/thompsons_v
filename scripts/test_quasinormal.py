from test import setup_script
setup_script(__file__)

from thompson.examples import *
from pprint import pprint

# import random
# seed = random.randint(0, 10**10)
# print(seed)
# random.seed(seed)
# psi, phi = random_conjugate_pair()

psi = aut_from_file('psi.aut')[0]
phi = aut_from_file('phi.aut')[0]

print('='*10)
print(psi.repeated_image('x a3 a3 a3', 4))
print(psi)
print(phi)
rho = psi.test_conjugate_to(phi)
print(rho)

# phi.quasinormal_basis())
# type_b, type_c = phi._split_basis()
# roots, graph = phi.equivalence_data(type_b, type_c)
# dump_graph(roots, graph)

# def plot_graph(G):
	# import matplotlib.pyplot as plt
	# edge_labels = { (u, v) : fmt_triple(data) for (u, v, data) in G.edges_iter(data=True) }
	# pos = nx.spring_layout(G)
	# nx.draw(G,
		# node_size = 3000,
		# pos = pos,
		# with_labels = True
	# )
	# nx.draw_networkx_edge_labels(G, pos, edge_labels)
	# plt.show()