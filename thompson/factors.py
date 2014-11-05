"""
.. testsetup::
	
	from pprint import pprint
	from collections import deque
	
	from thompson.factors import *
	from thompson.examples import *
"""
from collections import defaultdict, deque
from functools import partial
from io import StringIO
from itertools import chain, permutations
from pprint import pprint

import networkx as nx

from .automorphism import Automorphism
from .generators import Generators
from .word import Word

__all__ = ["PeriodicFactor", "InfiniteFactor"]

class AutomorphismFactor(Automorphism):
	"""An automorphism derived from a larger parent automorphism.
	
	:ivar domain_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	:ivar range_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	"""
	def __init__(self, domain, range, domain_relabeller=None, range_relabeller=None):
		"""In addition to creating an automorphism, we allow an optional relabelling homomorphism *relabeller* to be stored. This allows us to turn words in the factor back to words in the parent algebra.
		
		This class is returned by :meth:`~thompson.automorphism.Automorphism.free_factors` and doesn't need to be instantiated by the user.
		
		:raises LookupError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range*
		:raises TypeError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range* (as appropriate)
		
		.. seealso:: :meth:`Automorphism.__init__`, :meth:`~thompson.automorphism.Automorphism.free_factors`
		"""
		if domain_relabeller is None != range_relabeller is None:
			raise LookupError("Must either specify both relabellers or neither.")
		
		if domain_relabeller is not None and domain_relabeller.domain.signature != domain.signature:
			raise TypeError('Domain relabeller signature {} does not match automorphism\'s domain signature {}.'.format(
			  domain_relabeller.domain.signature, domain.signature))
		
		if range_relabeller is not None and range_relabeller.domain.signature != range.signature:
			raise TypeError('Range relabeller signature {} does not match automorphism\'s range signature {}.'.format(
			  range_relabeller.domain.signature, range.signature))
		
		super().__init__(domain, range)
		self.domain_relabeller = domain_relabeller
		self.range_relabeller = range_relabeller
	
	def relabel(self):
		"""Used to convert back to the parent algebra after doing computations in the derived algebra."""
		#TODO example
		if self.domain_relabeller is None or self.range_relabeller is None:
			raise AttributeError("This factor has not been assigned a relabeller.")
		return self.domain_relabeller.image_of_set(self.domain), self.range_relabeller.image_of_set(self.range)
	
	def __str__(self):
		output = StringIO()
		output.write(self._string_header())
		output.write("\nThis automorphism was derived from a parent automorphism.\n'x' and 'y' represent root words of the parent and current derived algebra, respectively.")
		
		domain_relabelled, range_relabelled = self.relabel()
		rows = self._format_table(
		    domain_relabelled, self.domain, self.range, range_relabelled,
		    sep = ['~>   ', '=>', '   ~>'], root_names = 'xyyx'
		)
		for row in rows:
			output.write('\n')
			output.write(row)
		return output.getvalue()

def get_factor_class(infinite):
	return InfiniteFactor if infinite else PeriodicFactor

class PeriodicFactor(AutomorphismFactor):
	r"""A purely periodic free factor which has been extracted from another component.
	
		>>> print(example_5_9_p)
		PeriodicFactor: V(2, 1) -> V(2, 1) specified by 7 generators (after expansion and reduction).
		This automorphism was derived from a parent automorphism.
		'x' and 'y' represent root words of the parent and current derived algebra, respectively.
		x1 a1 a1 a1 ~>    y1 a1 a1    => y1 a1 a2 a1    ~> x1 a1 a1 a2
		x1 a1 a1 a2 ~>    y1 a1 a2 a1 => y1 a1 a2 a2    ~> x1 a1 a2   
		x1 a1 a2    ~>    y1 a1 a2 a2 => y1 a1 a1       ~> x1 a1 a1 a1
		x1 a2 a1 a1 ~>    y1 a2 a1 a1 => y1 a2 a1 a2    ~> x1 a2 a1 a2
		x1 a2 a1 a2 ~>    y1 a2 a1 a2 => y1 a2 a1 a1    ~> x1 a2 a1 a1
		x1 a2 a2 a1 ~>    y1 a2 a2 a1 => y1 a2 a2 a2    ~> x1 a2 a2 a2
		x1 a2 a2 a2 ~>    y1 a2 a2 a2 => y1 a2 a2 a1    ~> x1 a2 a2 a1
		>>> sorted(example_5_9_p.cycle_type)
		[2, 3]
		>>> from pprint import pprint
		>>> #Two orbits of size 2, one orbit of size 3
		>>> pprint(example_5_9_p.multiplicity)
		{2: 2, 3: 1}
		
	
	:ivar multiplicity: a mapping :math:`d \mapsto m_\phi(d, X_\phi)` where :math:`\phi` is the current automorphism and :math:`X_\phi` is the :meth:`quasi-normal basis <thompson.automorphism.Automorphism.quasinormal_basis>` for :math:`\phi`.
	:ivar cycle_type: the set :math:`\{d \in \mathbb{N} : \text{$\exists$ an orbit of length $d$.}\}`
	"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		#See discussion before lemma 5.3.
		#Should have X_P = Y_P = Z_P
		assert self.quasinormal_basis() == self.domain
		assert self.quasinormal_basis().minimal_expansion_for(self) == self.domain
		
		#Need to have constructed the QNF basis here.
		self._setup_cycle_type()
		self._setup_multiplicities()
	
	def _setup_cycle_type(self):
		#see definition 5.8
		self.cycle_type = {type.data for type in self._qnf_orbit_types.values()}
	
	def _setup_multiplicities(self):
		counts = defaultdict(int)
		for orbit in self._qnf_orbit_types.values():
			orbit_size = orbit.data
			counts[orbit_size] += 1
		for orbit_size, occurances in counts.items():
			assert occurances % orbit_size == 0
			counts[orbit_size] = occurances // orbit_size
		self.multiplicity = dict(counts)
	
	def enumerate_orbits(self, basis):
		r"""Enumerates the periodic orbits of the current automorphism's quasinormal_basis. Returns a dictionary *orbits_by_size*. Each value ``orbits_by_size[d]`` is a list of the orbits of size *d*. Orbits themselves are represented as lists of :class:`Words <thompson.word.Word>`.
		
			>>> def display_orbits(orbits_by_size):
			... 	for key in sorted(orbits_by_size):
			... 		print('Orbits of length', key)
			... 		for list in orbits_by_size[key]:
			... 			print('...', *list, sep=' -> ', end=' -> ...\n')
			>>> orbits_by_size = example_5_9_p.enumerate_orbits(example_5_9_p.quasinormal_basis())
			>>> display_orbits(orbits_by_size)
			Orbits of length 2
			... -> x1 a2 a1 a1 -> x1 a2 a1 a2 -> ...
			... -> x1 a2 a2 a1 -> x1 a2 a2 a2 -> ...
			Orbits of length 3
			... -> x1 a1 a1 -> x1 a1 a2 a1 -> x1 a1 a2 a2 -> ...
		"""
		#TODO another doctest
		orbits_by_size = defaultdict(deque)
		already_seen = set()
		for gen in basis:
			if gen in already_seen:
				continue
			type, images = self.orbit_type(gen, basis)
			length = type.data
			images = [images[i] for i in range(length)]
			already_seen.update(images)
			orbits_by_size[length].append(images)
		return dict(orbits_by_size)
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely periodic automorphisms are periodic by examining their orbits.
		
			>>> psi_p = example_5_12_psi_p; phi_p = example_5_12_phi_p
			>>> rho_p = psi_p.test_conjugate_to(phi_p)
			>>> print(rho_p)
			AutomorphismFactor: V(2, 1) -> V(2, 1) specified by 6 generators (after expansion and reduction).
			This automorphism was derived from a parent automorphism.
			'x' and 'y' represent root words of the parent and current derived algebra, respectively.
			x1 a1 a1 a1 a1 ~>    y1 a1 a1    => y1 a1 a2       ~> x1 a1 a2   
			x1 a1 a1 a1 a2 ~>    y1 a1 a2    => y1 a2 a2       ~> x1 a2 a2   
			x1 a1 a1 a2    ~>    y1 a2 a1 a1 => y1 a1 a1 a1    ~> x1 a1 a1 a1
			x1 a1 a2       ~>    y1 a2 a1 a2 => y1 a2 a1 a1    ~> x1 a2 a1 a1
			x1 a2 a1       ~>    y1 a2 a2 a1 => y1 a1 a1 a2    ~> x1 a1 a1 a2
			x1 a2 a2       ~>    y1 a2 a2 a2 => y1 a2 a1 a2    ~> x1 a2 a1 a2
			>>> rho_p * phi_p == psi_p * rho_p
			True
		
		.. seealso:: This implements algorithm 5.13 of the paper---see section 5.3.
		"""
		# todo another doctest
		if not isinstance(other, PeriodicFactor):
			raise TypeError('Other automorphism must be a PeriodicFactor.')
		
		# 1. The quasi-normal bases are constructed in initialisation.
		# Note that the QNF basis should be just the domain. TODO Checkme
		
		# 2. Check that the cycle types are the same.
		if self.cycle_type != other.cycle_type:
			return None
		
		#3. Check that the multiplicites are congruent.
		modulus = self.signature.arity - 1
		for d in self.cycle_type:
			if self.multiplicity[d] % modulus != other.multiplicity[d] % modulus:
				return None
		
		# 4. Expand bases until the orbits multiplicites are the same
		s_orbits_of_size = self.enumerate_orbits(self.quasinormal_basis())
		o_orbits_of_size = other.enumerate_orbits(other.quasinormal_basis())

		for d in self.cycle_type:
			assert len(s_orbits_of_size[d]) == self.multiplicity[d]
			assert len(o_orbits_of_size[d]) == other.multiplicity[d]
			
			expansions_needed = (self.multiplicity[d] - other.multiplicity[d]) // modulus
			
			if expansions_needed > 0:
				other.expand_orbits(o_orbits_of_size[d], expansions_needed)
			elif expansions_needed < 0:
				self.expand_orbits(s_orbits_of_size[d], -expansions_needed)
			
			assert len(s_orbits_of_size[d]) == len(o_orbits_of_size[d])
		
		domain = Generators(self.signature)
		range  = Generators(self.signature)
		for d in self.cycle_type:
			for s_orbit, o_orbit in zip(s_orbits_of_size[d], o_orbits_of_size[d]):
				for s_word, o_word in zip(s_orbit, o_orbit):
					domain.append(s_word)
					range.append(o_word)
		
		rho = AutomorphismFactor(domain, range, self.domain_relabeller, other.range_relabeller)
		return rho
	
	@staticmethod
	def expand_orbits(deque, num_expansions):
		r"""Takes a *deque* whose elements are a list of words. The following process is repeated *num_expansions* times.
		
		1. Pop an orbit :math:`\mathcal{O} = [o_1, \dotsc, o_d]` from the left of the deque.
		2. For :math:`i = 1, \dotsc, n` (where :math:`n` is the arity):
			
			a. Compute the new orbit :math:`\mathcal{O_i} =[o_1\alpha_i, \dotsc, o_d\alpha_i]`
			b. Append this orbit to the right of *deque*.
		
		Once complete, the number of orbits in *deque* has increased by :math:`\text{(arity -1) $\times$ num_expansions}`
		"""
		for _ in range(num_expansions):
			orbit = deque.popleft()
			arity = orbit[0].signature.arity
			for i in range(1, arity+1):
				new_orbit = [w.alpha(i) for w in orbit]
				deque.append(new_orbit)

class InfiniteFactor(AutomorphismFactor):
	"""#todo docstring"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		#See discussion before lemma 5.3.
		#Domain should be Y_P
		assert self.quasinormal_basis().minimal_expansion_for(self) == self.domain
	
	def test_conjugate_to(self, other):
		#todo doctest and docstring
		if not isinstance(other, InfiniteFactor):
			raise TypeError('Other automorphism must be a InfiniteFactor.')
		
		#1. The QNF bases are computed automatically.
		#2. Compute the equivalence classes X_1, ... X_m of \equiv on self's QNF basis
		type_b, type_c = self._split_basis()
		print('** type B: **')
		pprint(type_b)
		print('** type C: **')
		pprint(type_c)
		roots, graph = self.equivalence_data(type_b, type_c)
		print('** Roots: **')
		pprint(roots)
		print('** Graph: **')
		dump_graph(roots, graph)
		#3. Find the initial and terminal elements of SI *other*-orbits.
		#4. Construct the sets R_i
		potential_endpoints = self.potential_image_endpoints(other, type_b)
		print('** Potential Endpoints: **')
		pprint(potential_endpoints)
		
		#5. Type B representitives for each class are stored in *roots*.
		#6. Iterate through all the automorhpisms which be conjugators.
		for rho in self.potential_conjugators(other, roots, graph, potential_endpoints, type_c):
			if self * rho == rho * other:
				return rho
		return None
	
	def _split_basis(self):
		"""Partition X into type B and type C parts."""
		#todo docstring and test
		type_b = {}
		type_c = {}
		basis = self.quasinormal_basis()
		for gen in basis:
			type, _ = self.orbit_type(gen, basis)
			if type.is_type('B'):
				type_b[gen] = type.data
			elif type.is_type('C'):
				type_c[gen] = type.data
			else:
				raise ValueError('Incorrect orbit type.')
		return type_b, type_c
	
	def equivalence_data(self, type_b, type_c):
		#todo doctest and string
		G = self._congruence_graph()
		self._simplify_graph(G, type_c)
		roots = self._reduce_to_forest(G)
		return roots, G
	
	def _congruence_graph(self):
		"""Form a graph whose vertex set is the QNF basis. The edges store the information which makes two vertices congruent."""
		basis = self.quasinormal_basis()
		min_exp = basis.minimal_expansion_for(self)
		endpts = self.semi_infinite_end_points()
		
		G = nx.DiGraph()
		orbit_generators = set(min_exp + endpts)
		
		#1. Add an edge for every direct conjugacy relationship.
		for gen in orbit_generators:
			type, images = self.orbit_type(gen, basis)
			for power, img in images.items():
				images[power] = basis.test_above(img)
			
			congruent_pairs = permutations(images.items(), 2)
			for (pow1, (head1, tail1)), (pow2, (head2, tail2)) in congruent_pairs:
				if head1 == head2:
					continue
				G.add_edge(head1, head2,
					  start_tail = tail1, power = pow2 - pow1, end_tail = tail2)
		return G
	
	@staticmethod
	def _simplify_graph(G, type_c):
		"""Removes all the type C words from this graph. Edges to and from type C words are removed, but the information they is still stored in the graph."""
		#2. Remove the type C elements.
		for word, type_b_data in type_c.items():
			replacement = type_b_data[1]
			replacement_in  = G[word][replacement]
			replacement_out = G[replacement][word]
			#Use the scheme of Lemma 5.24 to avoid type C words.
			for source, _, incoming in G.in_edges_iter(word, data=True):
				G.add_edge(source, replacement,
				  start_tail = incoming['start_tail'],
				  power      = incoming['power'] + replacement_in['power'],
				  end_tail   = replacement_in['end_tail'] + incoming['end_tail'])
				  
			for _, target, outgoing in G.out_edges_iter(word, data=True):
				G.add_edge(replacement, target,
				  start_tail = replacement_in['start_tail'] + incoming['start_tail'],
				  power      = replacement_out['power'] + incoming['power'],
				  end_tail   = outgoing['end_tail']) 
			
			G.remove_node(word)
		
		return [], G
	
	@staticmethod
	def _reduce_to_forest(G):
		"""Removes edges from G so that each connected component is a tree."""
		unseen_nodes = set(G.nodes_iter())
		unvisisted_edges = set(G.edges_iter())
		roots = []
		
		while unseen_nodes:
			root = unseen_nodes.pop()
			print('ROOT:', root)
			roots.append(root)
			for edge in nx.dfs_edges(G, root):
				unvisisted_edges.remove(edge)
				_, target = edge
				unseen_nodes.remove(target)
		
		to_remove = [e for e in G.edges_iter() if e in unvisisted_edges]
		G.remove_edges_from(to_remove)
		return roots
	
	def potential_image_endpoints(self, other, self_type_b):
		"""Let ``x`` be a type B word with respect to the current automorphism. This returns a mapping which takes ``x`` and produces the set of words ``w`` which are endpoints of *other*-orbits which have the same characteristic as ``x``.
		"""
		#todo doctest
		images_by_char = defaultdict(set)
		basis = other.quasinormal_basis()
		for word in other.semi_infinite_end_points():
			type, _ = other.orbit_type(word, basis)
			assert type.is_type('B')
			images_by_char[type.data].add(word)
		
		orbit_endpts = {word : images_by_char[char] for word, char in self_type_b.items()}
		return orbit_endpts
	
	def potential_conjugators(self, other, roots, graph, choices, type_c):
		#1. Flatten and glue the graph components together to form the ladder
		ladder = []
		for root in roots:
			ladder.extend(nx.dfs_preorder_nodes(graph, root))
		
		#2. Define the test function. See the docstring for maps_satisfying_choices
		deduce_images = partial(image_for_type_b, roots=roots, graph=graph, aut=other)
		
		#3. Now the hard part---the iteration.
		for images in maps_satisfying_choices(ladder, choices, deduce_images):
			domain = Generators(self.signature, sorted(images))
			range = Generators(self.signature, (images[d] for d in domain))
			
			#a. Add in type C images to domain and range
			for word, data in type_c.items():
				domain.append(word)
				img = image_for_type_c(word, data, images, other)
				range.append(img)
			
			#b. Does this define an Automorphism?
			try:
				pprint(domain)
				pprint(range)
				rho = AutomorphismFactor(domain, range, self.domain_relabeller, other.range_relabeller)
			except ValueError as e:
				print(e)
				continue
			yield rho

def maps_satisfying_choices(domain, choices, image_for):
	r"""Suppose we have a list of elements :math:`d` belonging to some list *domain*. We would like to efficiently enumerate the functions :math:`f\colon\text{domain} -> I` where :math:`I` is some set. We would also like these functions :math:`f` to abide by certain rules described below.
	
	For each such :math:`d` we have a set of *choices* to make: say :math:`\text{choices}(d) = \{c_{d,1}, \dots c_{d,n_d}\}`. Making a choice :math:`c_{d,i}` for :math:`d` will determine the set of potential images for :math:`d`. This set is allowed to depend on:
	
	- the element :math:`d` being mapped;
	- the choice :math:`c_{d, i}` which was made for d;
	- the images which have been proposed thus far for the elements preceeding :math:`d` in *domain*.
	
	We also ask that each choice produces at most one image for :math:`d`. This should be determined by the call ``image_for(d, choice, images)``, where images is the proposed mapping which has been constructed so far. If no image is available, this function should return ``None``.
	
	.. warning:: Make sure that ``None`` isn't a valid image for :math:`d`!
	
	:returns: yields mappings which satisfy all the constraints until it is no longer possible to do so.
	"""
	#2. Setup the state we need
	images = {}                             #mapping from ladder to images
	choice_iterators = [None] * len(domain) #iterators yielding the choices at each rung
	depth  = 0                              #current position on the ladder
	ascend = False                          #do we go up or down?
	
	while True:
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
				break
		
		#Select the next choice for this depth
		current = domain[depth]
		if choice_iterators[depth] is None:
			choice_iterators[depth] = iter(choices[current])
		try:
			choice = next(choice_iterators[depth])
		except GeneratorExit:
			ascend = True
			continue
		
		#Does this choice give us an image?
		print('trying the choice', current, '~>', choice)
		img = image_for(current, choice, images)
		if img is None:
			print('This doesn\'t work.')
			ascend = True
			continue
		print('WORKS, and', current, '->', img)
		images[current] = img
		ascend = False
		
		#If we've made it this far, we've chosen an image. Try to go down the ladder.
		if depth + 1 == len(domain):
			print("I YIELD CONTROL TO YOU")
			yield images
			ascend = True
		else:
			depth += 1

def image_for_type_b(word, chosen_endpoint, images, roots, graph, aut):
	"""Is it possible to map *word* into the *aut*-orbit ending with *chosen_endpoint* when we have already constructed a mapping described by *images*?"""
	if word in roots:
		return chosen_endpoint
	
	predecessor, _, edge_data = next(graph.in_edges_iter(word, data=True))
	predecessor_image = images[predecessor]
	
	u = chosen_endpoint.extend(edge_data['end_tail'])
	v = predecessor_image.extend(edge_data['start_tail'])
	solns = aut.share_orbit(u, v)
	if solns.is_empty():
		return None
	assert not solns.is_sequence()
	return aut.repeated_image(chosen_endpoint, solns.base + edge_data['power'])

def image_for_type_c(word, type_b_data, images, aut):
	power, gen, tail = type_b_data
	w = images[gen].extend(tail)
	return aut.repeated_image(w, -power)


###for debugging only
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




