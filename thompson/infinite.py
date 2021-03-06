"""
.. testsetup::
	
	from collections import deque
	from fractions   import gcd
	from pprint      import pprint
	
	from thompson          import *
	from thompson.word     import format
	from thompson.infinite import *
"""
from collections import defaultdict
from functools   import partial
from io          import StringIO
from itertools   import chain, permutations

from .number_theory import prod
from .word          import Word, root, format
from .generators    import Generators
from .homomorphism  import format_table
from .automorphism  import Automorphism
from .mixed         import MixedAut

__all__ = ["InfiniteAut"]

def get_factor_class(infinite):
	return InfiniteAut if infinite else PeriodicAut

class InfiniteAut(Automorphism):
	"""A purely infinite automorphism which may have been extracted from a mixed automorphism.
	
		>>> print(load_example('example_5_3').free_factors()[1])
		InfiniteAut: V(2, 1) -> V(2, 1) specified by 6 generators (after expansion and reduction).
		This automorphism was derived from a parent automorphism.
		'x' and 'y' represent root words of the parent and current derived algebra, respectively.
		x1 a1 a1 a1 a1 ~>    y1 a1 a1 a1 => y1 a1 a1       ~> x1 a1 a1 a1   
		x1 a1 a1 a1 a2 ~>    y1 a1 a1 a2 => y1 a1 a2 a1    ~> x1 a1 a1 a2 a1
		x1 a1 a1 a2    ~>    y1 a1 a2    => y1 a1 a2 a2    ~> x1 a1 a1 a2 a2
		x1 a2 a1       ~>    y1 a2 a1    => y1 a2 a1 a1    ~> x1 a2 a1 a1   
		x1 a2 a2 a1    ~>    y1 a2 a2 a1 => y1 a2 a1 a2    ~> x1 a2 a1 a2   
		x1 a2 a2 a2    ~>    y1 a2 a2 a2 => y1 a2 a2       ~> x1 a2 a2      
		>>> phi = random_infinite_automorphism()
		>>> phi.order
		inf
		>>> len(phi.cycle_type)
		0
		>>> len(phi.characteristics) != 0
		True
	"""
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely infinite automorphisms are conjugate by breaking down the quasi-normal basis into :meth:`equivalence classes <equivalence_data>`.
		
			>>> for psi_i, phi_i in (
			...   load_example_pair('example_5_26'),
			...   load_example_pair('inf_conj'),
			...   random_conjugate_infinite_pair()):
			... 	rho_i = psi_i.test_conjugate_to(phi_i)
			... 	print('Conjugate:', rho_i is not None, ' Conjugator works:', rho_i * phi_i == psi_i * rho_i)
			Conjugate: True  Conjugator works: True
			Conjugate: True  Conjugator works: True
			Conjugate: True  Conjugator works: True
			
			>>> f, g = (load_example('not_conjugate_' + c) for c in 'fg')
			>>> f.is_conjugate_to(g)
			False
		
		.. doctest::
			:hide:
			
			>>> rho = load_example('nathan_pond_example').test_conjugate_to(load_example('nathan1_example'))
			>>> rho is None
			True
			>>> f, g, h = (load_example('olga_' + c) for c in 'fgh')
			>>> rho = f.test_conjugate_to(g)
			>>> rho is not None
			True
			>>> f * rho == rho * g
			True
			>>> f * h == h * g
			True
			
			>>> d, e = ( load_example('two_flow_components_' + c) for c in 'de')
			>>> d.is_conjugate_to(e)
			True
		
		.. seealso:: This implements Algorithm :paperref:`alg:reginf` of the paper---see Section :paperref:`subsectionRI`.
		"""
		#todo another doctest.
		if not isinstance(other, InfiniteAut):
			return None
		#1. The QNF bases are computed automatically.
		#2. Compute the equivalence classes X_1, ... X_m of \equiv on self's QNF basis
		type_b, type_c = self._split_basis()
		roots, graph = self.equivalence_data(type_b, type_c)
		verify_graph(self, roots, graph)
		#3. Find the initial and terminal elements of SI *other*-orbits.
		#4. Construct the sets R_i
		potential_endpoints = other.potential_image_endpoints(type_b)
		
		#5. Type B representitives for each class are stored in *roots*.
		#6. Iterate through all the automorphisms which could be conjugators.
		for rho in self._potential_conjugators(other, roots, graph, potential_endpoints, type_c):
			if self * rho == rho * other:
				return rho
		return None
	
	def _split_basis(self):
		"""Partition X into type B and type C parts."""
		#todo docstring and test
		type_b = {}
		type_c = {}
		basis = self.quasinormal_basis
		for gen in basis:
			type, _, type_b_data = self.orbit_type(gen, basis)
			if type.is_type_B():
				type_b[gen] = type.characteristic
			elif type.is_type_C():
				type_c[gen] = type_b_data
			else:
				raise ValueError('Incorrect orbit type.')
		return type_b, type_c
	
	def equivalence_data(self, type_b, type_c):
		r"""Let :math:`X` be the quasi-normal basis for the current automorphism :math:`\psi`. We can define an equivalence relation :math:`\equiv` on :math:`X` by taking the non-transitive relation
		
		.. math:: x \equiv y \iff \exists k, \Gamma, \Delta : x\Gamma\psi^k = y\Delta
		
		and allowing it to generate an equivalence relation. This method returns a pair *(roots, graph)* where
		
		- *graph* is a :func:`~nx:networkx.DiGraph`
		
			- vertices are type B words in :math:`X`
			- directed edges :math:`x \to y` correspond to the direct relation :math:`x \equiv y`
			- the *graph* is a directed forest
		
		- *roots* is a list of type B words in :math:`X`
			
			- Each *root* is the root of a different tree in the forest *graph*.
		
		This information allows us to (attempt to) compute the images of :math:`X` under a conjugator given only the images of *roots*.
		
		.. seealso:: Definition :paperref:`EqCl` through to Lemma :paperref:`lem:effsplit`.
		"""
		G = self._congruence_graph()
		self._simplify_graph(G, type_c)
		roots = self._reduce_to_forest(G)
		return roots, G
	
	def _congruence_graph(self):
		"""Form a graph whose vertex set is the QNF basis. The edges store the information which makes two vertices congruent."""
		from networkx.classes.digraph import DiGraph
		basis = self.quasinormal_basis
		min_exp = basis.minimal_expansion_for(self)
		terminal, initial = self.semi_infinite_end_points()
		endpts = sorted(terminal + initial)
		G = DiGraph()
		orbit_generators = set(min_exp)
		orbit_generators.update(endpts)
		
		#1. Add an edge for every direct congruence relationship.
		for gen in orbit_generators:
			type, images, _ = self.orbit_type(gen, basis)
			for power, img in images.items():
				images[power] = basis.test_above(img)
			
			congruent_pairs = permutations(images.items(), 2)
			for (pow1, (head1, tail1)), (pow2, (head2, tail2)) in congruent_pairs:
				if head1 == head2:
					continue
				data = dict(start_tail = tail1, power = pow2 - pow1, end_tail = tail2)
				G.add_edge(head1, head2, data)
				assert self.repeated_image(head1.extend(tail1), pow2 - pow1) == head2.extend(tail2)
		return G
	
	# @staticmethod
	def _simplify_graph(self, G, type_c):
		"""Removes all the type C words from this graph. Edges to and from type C words are removed, but the information they is still stored in the graph."""
		#2. Remove the type C elements.
		for word, type_b_data in type_c.items():
			replacement = type_b_data['target']
			
			assert self.repeated_image(word, type_b_data['power']) == replacement.extend(type_b_data['end_tail'])
			
			#Use the scheme of Lemma 5.24 to avoid type C words.
			for source, _, incoming in G.in_edges_iter(word, data=True):
				if source == replacement:
					continue
				data = dict(
				  start_tail = incoming['start_tail'],
				  power      = incoming['power'] + type_b_data['power'],
				  end_tail   = type_b_data['end_tail'] + incoming['end_tail']
				)
				G.add_edge(source, replacement, data)
				assert self.repeated_image(source.extend(data['start_tail']), data['power']) == replacement.extend(data['end_tail'])
				  
			for _, target, outgoing in G.out_edges_iter(word, data=True):
				if target == replacement:
					continue
				data = dict(
				  start_tail = type_b_data['end_tail'] + outgoing['start_tail'],
				  power      = outgoing['power'] - type_b_data['power'],
				  end_tail   = outgoing['end_tail']
				)
				G.add_edge(replacement, target, data)
				assert self.repeated_image(replacement.extend(data['start_tail']), data['power']) == target.extend(data['end_tail'])
				
			G.remove_node(word)
		
		return [], G
	
	@staticmethod
	def _reduce_to_forest(G):
		"""Removes edges from G so that each connected component is a tree."""
		from networkx.algorithms.traversal.depth_first_search import dfs_edges
		unseen_nodes = set(G.nodes_iter())
		unvisisted_edges = set(G.edges_iter())
		roots = []
		
		while unseen_nodes:
			root = unseen_nodes.pop()
			roots.append(root)
			for edge in dfs_edges(G, root):
				unvisisted_edges.remove(edge)
				_, target = edge
				unseen_nodes.remove(target)
		
		to_remove = [e for e in G.edges_iter() if e in unvisisted_edges]
		G.remove_edges_from(to_remove)
		return roots
	
	def potential_image_endpoints(other, self_type_b):
		"""Let ``x`` be a type B word with respect to the current automorphism. This returns a mapping which takes ``x`` and produces the set of words ``w`` which are endpoints of *other*-orbits which have the same characteristic as ``x``.
		
		.. seealso:: The sets :math:`\mathcal R_i` of defintion :paperref:`def:potential_image_endpoints`.
		"""
		#todo doctest
		images_by_char = defaultdict(set)
		basis = other.quasinormal_basis
		terminal, initial = other.semi_infinite_end_points()
		for word in terminal + initial:
			type, _, _ = other.orbit_type(word, basis)
			if type.is_type_B():
				images_by_char[type.characteristic].add(word)
		
		orbit_endpts = {word : images_by_char[char] for word, char in self_type_b.items()}
		return orbit_endpts
	
	def _potential_conjugators(self, other, roots, graph, choices, type_c, verbose=False):
		log = get_logger(verbose)
		#1. Flatten and glue the graph components together to form the ladder
		from networkx.algorithms.traversal.depth_first_search import dfs_preorder_nodes
		ladder = []
		for root in roots:
			ladder.extend(dfs_preorder_nodes(graph, root))
		log("ladder is:")
		log(*ladder, sep="\n")
		
		#2. Define the test function. See the docstring for maps_satisfying_choices
		deduce_images = partial(image_for_type_b, roots=roots, graph=graph, aut=other)
		
		#3. Now the hard part---the iteration.
		for images in maps_satisfying_choices(ladder, choices, deduce_images, verbose):
			domain = Generators(self.signature, sorted(images))
			range = Generators(self.signature, (images[d] for d in domain))
			
			log("deducing images for type C words")
			#Add in type C images to domain and range
			for word, data in type_c.items():
				domain.append(word)
				img = image_for_type_c(word, data, images, other)
				range.append(img)
			
			try:
				rho = Automorphism(domain, range)
			except ValueError as e:
				log("This didn't define an automorphism")
				continue
			else:
				log("These images DID build an automorphism")
				rho.add_relabellers(self.domain_relabeller, other.range_relabeller)
				yield rho
	
	def find_power_conjugators(self, other, cheat=False):
		"""Tests two infinite factors to see if they are power conjugate. Yields minimal solutions :math:`(a, b, \rho)`.
		
			>>> example_6_8_psi, example_6_8_phi = load_example_pair('example_6_8')
			>>> solns = example_6_8_psi.find_power_conjugators(example_6_8_phi)
			>>> for a, b, rho in solns:
			... 	print(a, b)
			... 	assert example_6_8_psi ** a * rho == rho * example_6_8_phi ** b
			3 1
			-3 -1
		
		.. warning:: If the :meth:`~thompson.infinite.InfiniteAut.power_conjugacy_bounds` are reasonable large (say > 30), this method could potentially take a long time!
		
		.. seealso:: Algorithm :paperref:`alg:pcRI` of the paper.
		"""
		if not isinstance(other, InfiniteAut):
			raise StopIteration
		if cheat:
			from .examples.random import random_power_bounds
			bounds = random_power_bounds
		else:
			bounds = self.power_conjugacy_bounds(other)
		yield from self._test_power_conjugate_upto(other, *bounds, inverses=True)
	
	def test_power_conjugate_to(self, other, cheat=False):
		r"""Tests two infinite factors to see if they are power conjugate. If a solution exists, returns :math:`(a, b, \rho)` such that :math:`\rho^{-1}\psi^a\rho = \phi^b`. Otherwise returns None.
		
			>>> result = example_6_8_psi.test_power_conjugate_to(example_6_8_phi)
			>>> result is not None
			True
			>>> a, b, rho = result
			>>> (example_6_8_psi ** a) * rho == rho * (example_6_8_phi ** b)
			True
		
		.. warning:: If the :meth:`~thompson.infinite.InfiniteAut.power_conjugacy_bounds` are reasonable large (say > 30), this method could potentially take a long time!
		
		.. seealso:: Algorithm :paperref:`alg:pcRI` of the paper.
		"""
		try:
			return next(self.find_power_conjugators(other, cheat=cheat))
		except StopIteration:
			# print('No infinite conjugators')
			return None
	
	def power_conjugacy_bounds(self, other):
		"""Compute the bounds :math:`\hat a, \hat b` on powers which could solve the power conjugacy problem.
		
			>>> example_6_8_psi.power_conjugacy_bounds(example_6_8_phi)
			(9, 1)
		
		.. seealso:: Proposition :paperref:`divisors` and Corollary :paperref:`AJDCOROLLARYPC`.
		"""
		s_parts = self.minimal_partition()
		o_parts = other.minimal_partition()
		#If the root sets are different, there's no hope of power conjugacy.
		if s_parts.keys() != o_parts.keys():
			return 0, 0
		s_bound = self.compute_bounds(s_parts, o_parts)
		o_bound = self.compute_bounds(o_parts, s_parts)
		return s_bound, o_bound
	
	def minimal_partition(self):
		r"""Let :math:`\psi` be the current automorphism. This method partitions the characteristics :math:`M_\psi` into cells :math:`P_1 \sqcup \dots \sqcup P_L`, where
		- The multipliers :math:`\Gamma` all have the same :func:`~thompson.word.root`, for all :math:`(m, \Gamma)` in each :math:`P_i`.
		- :math:`L` is minimal with this property.
		
			>>> def print_partition(p):
			... 	for root in sorted(p.keys(), reverse=True):
			... 		print(format(root), end = ':')
			... 		for power, mult, _ in p[root]:
			... 			print(' ({}, {})'.format(power, format(mult)), end='')
			... 		print()
			>>> print_partition(example_6_8_psi.minimal_partition())
			a1: (-1, a1)
			a2: (1, a2)
			>>> print_partition(example_6_8_phi.minimal_partition())
			a1: (-1, a1 a1 a1)
			a2: (1, a2 a2 a2)
			>>> example_6_9_psi, example_6_9_phi = load_example_pair('example_6_9')
			>>> print_partition(example_6_9_phi.minimal_partition())
			a1: (-2, a1)
			a2: (1, a2)
			
		
		:returns: a dictionary of sets. The keys of this dictionary are the roots :math:`\sqrt\Gamma`; the values are the cells :math:`P_i`. An element of a cell looks like :math:`m, \Gamma, r` where :math:`m, \Gamma` is a characteristic and :math:`r` is the root power corresponding to :math:`\sqrt\Gamma`.
		
		.. seealso:: the discussion following Corollary :paperref:`AJDCOROLLARYPC`.
		"""
		chars = self.characteristics
		parts = defaultdict(set)
		for power, mult in chars:
			root_mult, root_power = root(mult)
			parts[root_mult].add((power, mult, root_power))
		return parts
	
	@staticmethod
	def compute_bounds(s_parts, o_parts):
		"""Computes the bounds :math:`\hat a, \hat b` (in terms of the partitions :math:`P, Q` given by *s_parts* and *o_parts*) as in eqns (14) and (15) of the paper.
		
			>>> def bounds(s, o):
			... 	P = s.minimal_partition()
			... 	Q = o.minimal_partition()
			... 	a_hat = InfiniteAut.compute_bounds(P, Q)
			... 	b_hat = InfiniteAut.compute_bounds(Q, P)
			... 	return a_hat, b_hat
			>>> bounds(example_6_8_psi, example_6_8_phi)
			(9, 1)
			>>> bounds(example_6_9_psi, example_6_9_phi)
			(1, 2)
		"""
		bound = 1
		for root in s_parts:
			P_i = s_parts[root]
			Q_i = o_parts[root]
			bound *= prod(abs(power) for power, mult, root_power in P_i) ** len(Q_i)
			bound *= prod(root_power for power, mult, root_power in Q_i) ** len(P_i)
		return bound

def maps_satisfying_choices(domain, choices, image_for, verbose=False):
	r"""Suppose we have a list of elements :math:`d` belonging to some list *domain*. We would like to efficiently enumerate the functions :math:`f\colon\text{domain} -> I` where :math:`I` is some set. We would also like these functions :math:`f` to abide by certain rules described below.
	
	For each such :math:`d` we have a set of *choices* to make: say :math:`\text{choices}(d) = \{c_{d,1}, \dots c_{d,n_d}\}`. Making a choice :math:`c_{d,i}` for :math:`d` will determine the set of potential images for :math:`d`. This set is allowed to depend on:
	
	- the element :math:`d` being mapped;
	- the choice :math:`c_{d, i}` which was made for d;
	- the images which have been proposed thus far for the elements preceeding :math:`d` in *domain*.
	
	We also ask that each choice produces at most one image for :math:`d`. This should be determined by the call ``image_for(d, choice, images)``, where images is the proposed mapping which has been constructed so far. If no image is available, this function should return ``None``.
	
	.. warning:: Make sure that ``None`` isn't a valid image for :math:`d`!
	
	:returns: yields mappings which satisfy all the constraints until it is no longer possible to do so.
	"""
	log = get_logger(verbose)
	#2. Setup the state we need
	images = {}                             #mapping from ladder to images
	choice_iterators = [None] * len(domain) #iterators yielding the choices at each rung
	depth  = 0                              #current position on the ladder
	ascend = False                          #do we go up or down?
	yield_images = False                    #have we found a possible mapping?
	
	while True:
		log("\ndepth", depth, "and going", "up" if ascend else "down")
		log("candidate images from ladder:")
		for key in domain:
			try:
				log(key, "->", images[key])
			except KeyError:
				break
		
		if yield_images:
			log("yielding these images")
			yield images
			yield_images = False
		
		#If the last choice we made didn't work (or we reached the bottom of the ladder) go up
		if ascend:
			current = domain[depth]
			try:
				del images[current]
			except KeyError:
				pass
			choice_iterators[depth] = None
			depth -= 1
			if depth < 0:
				log("all choices considered; ending loop")
				break
		ascend = False
		
		#Select the next choice for this depth
		current = domain[depth]
		log("where should we map", current)
		if choice_iterators[depth] is None:
			choice_iterators[depth] = iter(choices[current])
		try:
			choice = next(choice_iterators[depth])
			log("next choice:", choice) 
		except StopIteration:
			log("no more choices at depth", depth)
			ascend = True
			continue
		
		#Does this choice give us an image?
		img = image_for(current, choice, images)
		if img is None:
			log("no suitable image")
			continue
		if img in images.values():
			log("already mapped something to this image")
			continue
		log("corresponding image:", img)
		images[current] = img
		
		#If we've made it this far, we've chosen an image. Try to go down the ladder.
		if depth + 1 == len(domain):
			yield_images = True
			ascend = True
		else:
			depth += 1

def image_for_type_b(word, chosen_endpoint, images, roots, graph, aut):
	"""Is it possible to map *word* into the *aut*-orbit ending with *chosen_endpoint* when we have already constructed a mapping described by *images*?"""
	#In the paper's notation,
	#chosen_endpoint is w
	#word is x
	if word in roots:
		return chosen_endpoint
	
	predecessor, _, edge_data = next(graph.in_edges_iter(word, data=True))
	predecessor_image = images[predecessor]
	
	u = chosen_endpoint.extend(edge_data['end_tail'])     #w   Delta
	v = predecessor_image.extend(edge_data['start_tail']) #y_i Gamma
	solns = aut.share_orbit(u, v)
	if solns.is_empty():
		return None
	assert not solns.is_sequence(), (u, v, solns)
	return aut.repeated_image(chosen_endpoint, solns.base + edge_data['power'])


def image_for_type_c(word, type_b_data, images, aut):
	"""The end of Lemma 5.24 shows how we can calculate the image of a type C word given the images of all type B words."""
	power, gen, tail = type_b_data['power'], type_b_data['target'], type_b_data['end_tail']
	w = images[gen].extend(tail)
	return aut.repeated_image(w, -power)

#For debugging only
def fmt_triple(edge_data):
	return "({}, {}, {})".format(
		format(edge_data['start_tail']), edge_data['power'], format(edge_data['end_tail'])
	)

def print_edge(source, target, data):
	print("{} --> {} with data {}".format(source, target, fmt_triple(data)))

def dump_graph(roots, graph):
	from networkx.algorithms.traversal.depth_first_search import dfs_preorder_nodes
	for i, root in enumerate(roots):
		print('component', i, 'with root', root)
		print('type B elements:')
		for node in dfs_preorder_nodes(graph, root):
			print('\tEdges out of', node)
			for source, target in graph.out_edges_iter(node):
				data = graph[source][target]
				print('\t\tto', target, 'with data\n\t\t\t', fmt_triple(data))

def verify_graph(aut, roots, graph):
	from networkx.algorithms.traversal.depth_first_search import dfs_preorder_nodes
	for i, root in enumerate(roots):
		for node in dfs_preorder_nodes(graph, root):
			for source, target in graph.out_edges_iter(node):
				data = graph[source][target]
				assert aut.repeated_image(source.extend(data['start_tail']), data['power']) == target.extend(data['end_tail'])

def get_logger(verbose):
	if verbose:
		log = print
	else:
		def log(*args, **kwargs): pass
	return log