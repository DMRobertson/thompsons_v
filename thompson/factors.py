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

from .generators import Generators
from .homomorphism import format_table
from .automorphism import Automorphism
from .word import Word

__all__ = ["PeriodicFactor", "InfiniteFactor"]

class AutomorphismFactor(Automorphism):
	"""An automorphism derived from a larger parent automorphism.
	
	:ivar domain_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	:ivar range_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	"""
	#todo docstring for reduce
	def __init__(self, domain, range, domain_relabeller=None, range_relabeller=None, reduce=True):
		"""In addition to creating an automorphism, we allow optional relabelling homomorphisms to be stored. This allows us to turn words in the factor back into words in the parent algebra.
		
		Instances of this class are returned by :meth:`~thompson.automorphism.Automorphism.free_factor` and don't need to be instantiated by the user.
		
		:raises LookupError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range*
		:raises TypeError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range* (as appropriate)
		
		.. seealso:: :meth:`thompson.automorphism.Automorphism.__init__`, :meth:`~thompson.automorphism.Automorphism.free_factor`
		"""
		if domain_relabeller is None != range_relabeller is None:
			raise LookupError("Must either specify both relabellers or neither.")
		
		if domain_relabeller is not None and domain_relabeller.domain.signature != domain.signature:
			raise TypeError('Domain relabeller signature {} does not match automorphism\'s domain signature {}.'.format(
			  domain_relabeller.domain.signature, domain.signature))
		
		if range_relabeller is not None and range_relabeller.domain.signature != range.signature:
			raise TypeError('Range relabeller signature {} does not match automorphism\'s range signature {}.'.format(
			  range_relabeller.domain.signature, range.signature))
		
		super().__init__(domain, range, reduce)
		self.domain_relabeller = domain_relabeller
		self.range_relabeller = range_relabeller
	
	def relabel(self):
		r"""Used to convert back to the parent algebra after doing computations in the derived algebra.
		
		In the following example :meth:`~thompson.automorphism.Automorphism.test_conjugate_to` takes a pure periodic automorphism and extracts factors. A conjugator :math:`\rho` is produced by :meth:`the overridden version of this method <thompson.factors.PeriodicFactor.test_conjugate_to>`. Finally :math:`\rho` is relabelled back to the parent algebra.
		
		:raises AttributeError: if the factor has not been assigned any relabellers.
		
			>>> psi = example_5_12_psi; phi = example_5_12_phi
			>>> rho = psi.test_conjugate_to(phi)
			>>> print(rho)
			Automorphism: V(2, 1) -> V(2, 1) specified by 6 generators (after expansion and reduction).
			x1 a1 a1 a1 a1 -> x1 a1 a2   
			x1 a1 a1 a1 a2 -> x1 a2 a2   
			x1 a1 a1 a2    -> x1 a1 a1 a1
			x1 a1 a2       -> x1 a2 a1 a1
			x1 a2 a1       -> x1 a1 a1 a2
			x1 a2 a2       -> x1 a2 a1 a2
			>>> rho * phi == psi * rho
			True
		"""
		if self.domain_relabeller is None or self.range_relabeller is None:
			raise AttributeError("This factor has not been assigned relabellers.")
		return self.domain_relabeller.image_of_set(self.domain), self.range_relabeller.image_of_set(self.range)
	
	def __str__(self):
		if self.domain_relabeller is None:
			return super().__str__()
		output = StringIO()
		output.write(self._string_header())
		output.write("\nThis automorphism was derived from a parent automorphism.\n'x' and 'y' represent root words of the parent and current derived algebra, respectively.")
		
		domain_relabelled, range_relabelled = self.relabel()
		rows = format_table(
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
	r"""A purely periodic free factor which has been extracted from another automorphism.
	
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
		#See discussion before lemma 5.3: should have X_P = Y_P = Z_P
		assert self.quasinormal_basis() == self.domain
		assert self.quasinormal_basis().minimal_expansion_for(self) == self.domain
		
		self._setup_multiplicities()
	
	def _setup_multiplicities(self):
		#see definition 5.8
		self.cycle_type = set()
		counts = defaultdict(int)
		basis = self.quasinormal_basis()
		
		for gen in basis:
			type, _, _ = self.orbit_type(gen, basis)
			orbit_size = type.characteristic[0]
			self.cycle_type.add(orbit_size)
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
			type, images, _ = self.orbit_type(gen, basis)
			length = type.characteristic[0]
			images = [images[i] for i in range(length)]
			already_seen.update(images)
			orbits_by_size[length].append(images)
		return dict(orbits_by_size)
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely periodic automorphisms are conjugate by examining their orbits.
		
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
			
			>>> psi, phi = random_conjugate_pair()
			>>> psi_p = psi.free_factors()[0]
			>>> phi_p = phi.free_factors()[0]
			>>> rho_p = psi_p.test_conjugate_to(phi_p)
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
				expand_orbits(o_orbits_of_size[d], expansions_needed)
			elif expansions_needed < 0:
				expand_orbits(s_orbits_of_size[d], -expansions_needed)
			
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
	"""A purely infinite free factor which has been extracted from another automorphism.
	
		>>> print(example_5_3_i)
		InfiniteFactor: V(2, 1) -> V(2, 1) specified by 6 generators (after expansion and reduction).
		This automorphism was derived from a parent automorphism.
		'x' and 'y' represent root words of the parent and current derived algebra, respectively.
		x1 a1 a1 a1 a1 ~>    y1 a1 a1 a1 => y1 a1 a1       ~> x1 a1 a1 a1   
		x1 a1 a1 a1 a2 ~>    y1 a1 a1 a2 => y1 a1 a2 a1    ~> x1 a1 a1 a2 a1
		x1 a1 a1 a2    ~>    y1 a1 a2    => y1 a1 a2 a2    ~> x1 a1 a1 a2 a2
		x1 a2 a1       ~>    y1 a2 a1    => y1 a2 a1 a1    ~> x1 a2 a1 a1   
		x1 a2 a2 a1    ~>    y1 a2 a2 a1 => y1 a2 a1 a2    ~> x1 a2 a1 a2   
		x1 a2 a2 a2    ~>    y1 a2 a2 a2 => y1 a2 a2       ~> x1 a2 a2      
	"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		#See discussion before lemma 5.3.
		#Domain should be Y_P
		assert self.quasinormal_basis().minimal_expansion_for(self) == self.domain
	
	def test_conjugate_to(self, other):
		"""We can determine if two purely infinite automorphisms are conjugate by breaking down the :meth:`quasi-normal basis <thompson.automorphism.Automorphism.quasinormal_basis>` into :meth:`equivalence classes <equivalence_data>`.
		
		.. todo:: For the purposes of these doctests, need to have some way of making this output reproduceable. Will depend on how the equivalence graph is formed and iterated over.
		
			>>> psi_i = example_5_26_psi_i; phi_i = example_5_26_phi_i
			>>> rho_i = psi_i.test_conjugate_to(phi_i)
			>>> rho_i is not None
			True
			>>> rho_i * phi_i == psi_i * rho_i
			True
		
		.. seealso:: This implements algorithm 5.27 of the paper---see section 5.4.
		"""
		#todo another doctest.
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
		verify_graph(self, roots, graph)
		#3. Find the initial and terminal elements of SI *other*-orbits.
		#4. Construct the sets R_i
		print('** Potential Endpoints: **')
		potential_endpoints = other.potential_image_endpoints(type_b)
		pprint(potential_endpoints)
		
		#5. Type B representitives for each class are stored in *roots*.
		#6. Iterate through all the automorhpisms which be conjugators.
		for rho in self._potential_conjugators(other, roots, graph, potential_endpoints, type_c):
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
			type, _, type_b_data = self.orbit_type(gen, basis)
			if type.is_type_B():
				type_b[gen] = type.characteristic
			elif type.is_type_C():
				type_c[gen] = type_b_data
			else:
				raise ValueError('Incorrect orbit type.')
		return type_b, type_c
	
	def equivalence_data(self, type_b, type_c):
		r"""Let :math:`X` be the quasi-normal basis for the current automorphism :math:`\psi`. We can define an equivalence relation :math:`\equiv` on :math:`X` by taking the relation
		
		.. math:: x \equiv y \iff \exists k, \Gamma, \Delta : x\Gamma\psi^k = y\Delta
		
		and allowing it to generate an equivalence relation. This method returns a pair *(roots, graph)* where
		
		- *graph* is a `DiGraph <https://networkx.github.io/documentation/latest/reference/classes.digraph.html>`_
		
			- vertices are type B words in :math:`X`
			- directed edges :math:`x \to y` correspond to the relation :math:`x \equiv y`
			- the *graph* is a directed forest
		
		- *roots* is a list of type B words in :math:`X`
			
			- *roots* are the roots of the trees in *graph*
		
		This information allows us to (attempt to) compute the images of :math:`X` under a conjugator given only the images of *roots*.
		
		.. seealso:: Definition 5.17 through to Lemma 5.20; also section 2.2 of the paper.
		"""
		G = self._congruence_graph()
		self._simplify_graph(G, type_c)
		roots = self._reduce_to_forest(G)
		return roots, G
	
	def _congruence_graph(self):
		"""Form a graph whose vertex set is the QNF basis. The edges store the information which makes two vertices congruent."""
		basis = self.quasinormal_basis()
		min_exp = basis.minimal_expansion_for(self)
		terminal, initial = self.semi_infinite_end_points()
		endpts = sorted(terminal + initial)
		G = nx.DiGraph()
		orbit_generators = set(min_exp + endpts)
		
		#1. Add an edge for every direct conjugacy relationship.
		print('** edges **')
		for gen in orbit_generators:
			type, images, _ = self.orbit_type(gen, basis)
			for power, img in images.items():
				images[power] = basis.test_above(img)
			
			congruent_pairs = permutations(images.items(), 2)
			for (pow1, (head1, tail1)), (pow2, (head2, tail2)) in congruent_pairs:
				if head1 == head2:
					continue
				G.add_edge(head1, head2,
					  start_tail = tail1, power = pow2 - pow1, end_tail = tail2)
				assert self.repeated_image(head1.extend(tail1), pow2 - pow1) == head2.extend(tail2)
		return G
	
	# @staticmethod
	def _simplify_graph(self, G, type_c):
		"""Removes all the type C words from this graph. Edges to and from type C words are removed, but the information they is still stored in the graph."""
		print('** Simplify graph **')
		#2. Remove the type C elements.
		for word, type_b_data in type_c.items():
			replacement = type_b_data[1]
			replacement_in  = G[word][replacement]
			replacement_out = G[replacement][word]
			
			print('Type C: {}, Type B replacement: {}'.format(word, replacement, replacement_in))
			print_edge(word, replacement, replacement_in)
			assert self.repeated_image(word.extend(replacement_in['start_tail']), replacement_in['power']) == replacement.extend(replacement_in['end_tail'])
			print_edge(replacement, word, replacement_out)
			assert self.repeated_image(replacement.extend(replacement_out['start_tail']), replacement_out['power']) == word.extend(replacement_out['end_tail'])
			
			#Use the scheme of Lemma 5.24 to avoid type C words.
			for source, _, incoming in G.in_edges_iter(word, data=True):
				if source == replacement:
					continue
				data = dict(
				  start_tail = incoming['start_tail'],
				  power      = incoming['power'] + replacement_in['power'],
				  end_tail   = replacement_in['end_tail'] + incoming['end_tail']
				)
				G.add_edge(source, replacement, data)
				print('adding', end=' ');
				print_edge(source, replacement, data)
				assert self.repeated_image(source.extend(data['start_tail']), data['power']) == replacement.extend(data['end_tail'])
				  
			for _, target, outgoing in G.out_edges_iter(word, data=True):
				if target == replacement:
					continue
				data = dict(
				  start_tail = replacement_out['start_tail'] + outgoing['start_tail'],
				  power      = replacement_out['power'] + outgoing['power'],
				  end_tail   = outgoing['end_tail']
				)
				G.add_edge(replacement, target, data) 
				print('adding', end=' ');
				print_edge(replacement, target, data)
				assert self.repeated_image(replacement.extend(data['start_tail']), data['power']) == target.extend(data['end_tail'])
				
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
	
	def potential_image_endpoints(other, self_type_b):
		"""Let ``x`` be a type B word with respect to the current automorphism. This returns a mapping which takes ``x`` and produces the set of words ``w`` which are endpoints of *other*-orbits which have the same characteristic as ``x``.
		
		.. seealso:: The sets :math:`\mathcal R_i` of defintion 5.23.
		"""
		#todo doctest
		images_by_char = defaultdict(set)
		basis = other.quasinormal_basis()
		terminal, initial = other.semi_infinite_end_points()
		for word in terminal + initial:
			type, _, _ = other.orbit_type(word, basis)
			if type.is_type_B():
				print('*', word, type)
				images_by_char[type.characteristic].add(word)
			else:
				pass
				print(' ', word, type)
		
		orbit_endpts = {word : images_by_char[char] for word, char in self_type_b.items()}
		return orbit_endpts
	
	def _potential_conjugators(self, other, roots, graph, choices, type_c):
		#1. Flatten and glue the graph components together to form the ladder
		ladder = []
		for root in roots:
			ladder.extend(nx.dfs_preorder_nodes(graph, root))
		
		#2. Define the test function. See the docstring for maps_satisfying_choices
		deduce_images = partial(image_for_type_b, roots=roots, graph=graph, aut=other)
		
		#3. Now the hard part---the iteration.
		for images in maps_satisfying_choices(ladder, choices, deduce_images):
			try:
				domain = Generators(self.signature, sorted(images))
				range = Generators(self.signature, (images[d] for d in domain))
				
				#Add in type C images to domain and range
				for word, data in type_c.items():
					domain.append(word)
					img = image_for_type_c(word, data, images, other)
					range.append(img)
				rho = AutomorphismFactor(domain, range, self.domain_relabeller, other.range_relabeller)
			except ValueError:
				#We haven't specified an automorphism.
				continue
			else:
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
	print(' ** Maps satisfying choices ** ')
	#2. Setup the state we need
	images = {}                             #mapping from ladder to images
	choice_iterators = [None] * len(domain) #iterators yielding the choices at each rung
	depth  = 0                              #current position on the ladder
	ascend = False                          #do we go up or down?
	print(domain)
	while True:
		#If the last choice we made didn't work (or we reached the bottom of the ladder) go up
		if ascend:
			print('Going up')
			current = domain[depth]
			try:
				del images[current]
			except KeyError:
				pass
			choice_iterators[depth] = None
			depth -= 1
			if depth < 0:
				print('break')
				break
		ascend = False
		
		#Select the next choice for this depth
		current = domain[depth]
		if choice_iterators[depth] is None:
			print('create choice iterator for depth', depth)
			choice_iterators[depth] = iter(choices[current])
		try:
			choice = next(choice_iterators[depth])
			print('trying the choice', current, '~> orbit ending in', choice)
		except StopIteration:
			print('No more choices available')
			ascend = True
			continue
		
		#Does this choice give us an image?
		img = image_for(current, choice, images)
		if img is None:
			print('This doesn\'t work.')
			continue
		print('WORKS, and', current, '->', img)
		images[current] = img
		
		#If we've made it this far, we've chosen an image. Try to go down the ladder.
		if depth + 1 == len(domain):
			print("I YIELD CONTROL TO YOU")
			yield images
			ascend = True
		else:
			print('going down')
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
	power, gen, tail = type_b_data
	w = images[gen].extend(tail)
	return aut.repeated_image(w, -power)

###for debugging only
def fmt_triple(edge_data):
	return "({}, {}, {})".format(
		format(edge_data['start_tail']), edge_data['power'], format(edge_data['end_tail'])
	)

def print_edge(source, target, data):
	print("{} --> {} with data {}".format(source, target, fmt_triple(data)))

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

def verify_graph(aut, roots, graph):
	import networkx as nx
	for i, root in enumerate(roots):
		for node in nx.dfs_preorder_nodes(graph, root):
			for source, target in graph.out_edges_iter(node):
				data = graph[source][target]
				assert aut.repeated_image(source.extend(data['start_tail']), data['power']) == target.extend(data['end_tail'])





