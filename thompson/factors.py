"""
.. testsetup::
	
	from thompson.factors import *
	from thompson.examples import *
"""
from collections import defaultdict, deque
from io import StringIO
from itertools import chain, permutations

import networkx as nx

from .automorphism import Automorphism
from .generators import Generators
from .word import Word

class AutomorphismFactor(Automorphism):
	"""An automorphism derived from a larger parent automorphism.
	
	:ivar domain_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	:ivar range_relabeller: the homomorphism which maps relabelled words back into the original algebra that this automorphism came from.
	"""
	def __init__(self, domain, range, domain_relabeller=None, range_relabeller=None):
		"""In addition to creating an automorphism, we allow an optional relabelling homomorphism *relabeller* to be stored. This allows us to turn words in the factor back to words in the parent algebra.
		
		This class is returned by :meth:`~thompson.automorphism.Automorphism.free_factors` and doesn't need to be instantiated by the user.
		
		:raises ValueError: if a relabeller is provided and the relabeller's signature does not match that of *domain* or *range*
		
		.. seealso:: :meth:`Automorphism.__init__`, :meth:`~thompson.automorphism.Automorphism.free_factors`
		"""
		if domain_relabeller is None != range_relabeller is None:
			raise ValueError("Must either specify both relabellers or neither.")
		
		if domain_relabeller is not None and domain_relabeller.domain.signature != domain.signature:
			raise ValueError('Domain relabeller signature {} does not match automorphism\'s domain signature {}.'.format(
			  domain_relabeller.domain.signature, domain.signature))
		
		if range_relabeller is not None and range_relabeller.domain.signature != range.signature:
			raise ValueError('Range relabeller signature {} does not match automorphism\'s range signature {}.'.format(
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

class PeriodicFactor(AutomorphismFactor):
	r"""A purely periodic free factor which has been extracted from another component.
	
		>>> print(example_5_9_p)
		PeriodicFactor: V(2, 1) -> V(2, 1) specified by 7 generators (after reduction).
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
			AutomorphismFactor: V(2, 1) -> V(2, 1) specified by 6 generators (after reduction).
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
		#1. The QNF bases are computed automatically.
		#2. Compute the equivalence classes X_1, ... X_m of \equiv on self's QNF basis
		classes = self.equivalence_classes()

		
		#3. Find the initial and terminal elements of SI *other*-orbits.
		endpts = other.semi_infinite_end_points()
		print(Generators.__str__(endpts))
		
		#4. Construct the sets R_i
		pprint(
			self.type_b_images(other, classes, endpts))
	
	def equivalence_classes(self):
		"""Partitions the QNF basis into equivalence classes. A class is represented as a triple *(graph, root, type_c_data)*, where:
		- the graph's vertices are the type B words in the class.
		- the graph's edges stores how the type B words are related
		- root is a chosen vertex which preceeds all others in the graph
		- type_c_data is a dictionary. Keys are the type C elements of this class. Values are the type B data recorded by the OrbitType class.
		"""
		type_b, type_c = self._split_basis()
		G = self._congruence_graph(type_b)
		components, roots = self._tree_components(G)
		#Before returning, we need to reinclude the type C elements!
		type_c_data = []
		for comp in components:
			restriction = {head : data for head, data in type_c.items() if data[1] in comp}
			type_c_data.append(restriction)
		return list(zip(components, roots, type_c_data))
	
	def _split_basis(self):
		"""Partition X into type B and type C parts."""
		#todo docstring and test
		type_b = []
		type_c = {}
		basis = self.quasinormal_basis()
		for gen in basis:
			type, _ = self.orbit_type(gen, basis)
			if type.is_type('B'):
				type_b.append(gen)
			elif type.is_type('C'):
				type_c[gen] = type.data
			else:
				raise ValueError('Incorrect orbit type.')
		return type_b, type_c
	
	def _congruence_graph(self, type_b):
		"""Form the graph whose vertices are type B elements of the QNF basis and edges store the information which makes two vertices congruent."""
		basis = self.quasinormal_basis()
		min_exp = basis.minimal_expansion_for(self)
		endpts = self.semi_infinite_end_points()
		
		G = nx.DiGraph()
		orbit_generators = set(min_exp + endpts)
		
		for gen in orbit_generators:
			type, images = self.orbit_type(gen, basis)
			type_b_images = {}
			for power, image in images.items():
				head, tail = basis.test_above(image)
				if head in type_b:
					type_b_images[power] = (head, tail)
			
			congruent_pairs = permutations(type_b_images.items(), 2)
			for (pow1, (head1, tail1)), (pow2, (head2, tail2)) in congruent_pairs:
				if head1 == head2:
					continue #no loops in this graph, please!
				# print("[{}] {} PSI^{} = [{}] {}".format(
				  # head1, Word.__str__(tail1), pow2-pow1, head2, Word.__str__(tail2)))
				G.add_edge(head1, head2,
					  start_tail = tail1, power = pow2 - pow1, end_tail = tail2)
		return G
		
	def _tree_components(self, G):
		"""Breaks down the _congruence_graph() into connected components, then removes extra stuff to get """
		unseen = set(G.nodes_iter())
		components = []
		roots = []
		
		while unseen:
			root = unseen.pop()
			T = nx.dfs_tree(G, root)
			components.append(T)
			roots.append(root)
			#Copy over the edge data
			for (u, v) in T.edges_iter():
				T[u][v] = G[u][v]
		return components, roots
	
	def semi_infinite_end_points(self):
		#todo docstring and test
		basis = self.quasinormal_basis()
		min_expansion = basis.minimal_expansion_for(self)
		img_expansion = self.image_of_set(min_expansion)
		terminal = basis.descendants_above(min_expansion)
		initial  = basis.descendants_above(img_expansion)
		return initial + terminal
	
	'''def type_b_images(self, other, classes, endpts):
		#1. Put each endpoint into a subalgebra V_i
		endpts_in_subalgebra = [set() for cls in classes]
		words_above_class = [set(cls.simple_words_above()) for cls in classes]
		for endpt in endpts:
			for i, cls in enumerate(classes):
				#1. Is the endpoint below any of the classes?
				if cls.is_above(endpt) or endpt in words_above_class[i]:
					endpts_in_subalgebra[i].add(endpt)
					break
			raise ValueError('Could not put {} into a subalgebra.'.format(
			  endpt))
		
		sbasis = self.quasinormal_basis()
		obasis = other.quasinormal_basis()
		type_b_images = []
		
		#2. For each subalgebra:
		for i, cls in enumerate(classes):
			#a. Map characteristics to the set of endpts in the subalgebra with those characteristics
			endpts_by_char = {}
			for endpt in endpts_in_subalgebra[i]:
				type, _ = other.orbit_type(endpt, obasis)
				try:
					endpts_by_char[type.data].add(endpt)
				except KeyError:
					endpts_by_char[type.data] = set(endpt)
			
			#b. Map type B words to the set of endpts with the same characteristic
			images = {}
			for word in cls:
				type, _ = self.orbit_type(word, sbasis)
				if not type.is_type('B'):
					continue
				images[word] = endpts_by_char[type.data]
			type_b_images.append(images)
		return type_b_images'''

def get_factor_class(infinite):
	return InfiniteFactor if infinite else PeriodicFactor