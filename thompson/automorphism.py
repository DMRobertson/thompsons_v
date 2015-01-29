"""
.. testsetup::
	
	from thompson.number_theory import gcd
	from thompson.word          import from_string
	from thompson.orbits        import print_component_types
	from thompson.automorphism  import *
	from thompson.examples      import *
"""

from copy      import copy
from itertools import chain, product

from .number_theory import lcm
from .word          import Word, free_monoid_on_alphas, format
from .generators    import Generators
from .homomorphism  import Homomorphism
from .orbits        import ComponentType, SolutionSet

___all__ = ["Automorphism"]

class Automorphism(Homomorphism):
	r"""
	Generic attributes:
	:ivar signature: The :class:`~thompson.word.Signature` shared by the generating sets domain and range.
	:ivar quasinormal_basis: See :meth:`compute_quasinormal_basis`.
	:ivar pond_banks: A list of tuples :math:`(\ell, k, r)` such that :math:`(\ell, r)` are banks of a pond with :math:`\ell\phi^k = r`.
	
	Periodic attributes:
	:ivar multiplicity: a mapping :math:`d \mapsto m_\phi(d, X_\phi)` where :math:`\phi` is the current automorphism and :math:`X_\phi` is the :meth:`quasi-normal basis <thompson.mixed.MixedAut.quasinormal_basis>` for :math:`\phi`.
	:ivar cycle_type: the set :math:`\{d \in \mathbb{N} : \text{$\exists$ an orbit of length $d$.}\}`
	:ivar order: the smallest positive number :math:`n` for which :math:`phi^n` is the identity. (This is the lcm of the cycle type.) If no such :math:`n` exists, the order is :math:`\infty`.
	
	Infinite attributes:
	:ivar characteristics: the set of characteristics :math:`(m, \Gamma)` of this automorphism.
	
	"""
	#Initialisation
	def __init__(self, domain, range, reduce=True):
		r"""Creates an automorphism mapping the given *domain* basis to the *range* basis in the given order.
		
		.. math:: \text{domain}_{\,i} \mapsto \text{range}_{\,i}
		
		:raises TypeError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis :meth:`isn't actually a basis <thompson.generators.Generators.is_basis>`.
		
		.. seealso:: :meth:`thompson.homomorphism.Homomorphism.__init__`
		"""
		#Check to see that the domain and range algebras are the same
		if domain.signature != range.signature:
			raise TypeError('Domain {} and range {} have different signatures.'.format(
			  domain.signature, range.signature))
		
		#Check to see that range is a basis for the range algebra
		i, j = range.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		self._inv = {}
		super().__init__(domain, range, reduce)
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{}. Missing elements are {}.".format(
			  range.signature, [format(x) for x in missing]))
		
		#Everthing seems to be okay here.
		self.signature = domain.signature
		
		#Setup the inverse map
		for root in Generators.standard_basis(domain.signature):
			self._image_simple_above_domain(root, self.range.signature, self.domain.signature, self._inv)
		
		self.compute_quasinormal_basis()
	
	#Computing images
	def _set_image(self, d, r, sig_in, sig_out, cache):
		if cache is not self._map and cache is not self._inv:
			raise ValueError("Incorrect cache provided.")
		super()._set_image(d, r, sig_in, sig_out, cache)
		cache = self._inv if cache is self._map else self._map
		cache[r] = d
	
	def image(self, key, inverse=False):
		"""If *inverse* is True, the inverse of the current automorphism is used to map *key* instead. Otherwise this method delegates to :meth:`Homomorphism.image <thompson.homomorphism.Homomorphism.image>`.
		
		Examples of finding inverse images:
		
			>>> print(example_4_25.image('x1 a2 a2', inverse=True))
			x1 a2 a2 a1 a1
			>>> print(example_4_25.image('x1 a1 a1 a2 a2 a1', inverse=True))
			x1 a2 a1 a2 a1
			>>> print(example_4_25.image('x a2', inverse=True))
			x1 a2 a2 a2 x1 a2 a2 a1 a1 L
			>>> print(example_4_25.image('x a2 a2 x a1 a2 L', inverse=True))
			x1 a2 a2 a1
		"""
		if inverse:
			return super().image(key, self.range.signature, self.domain.signature, self._inv)
		else:
			return super().image(key)
	
	def image_of_set(self, set, inverse=False):
		"""If *inverse* is True, the inverse of the current automorphism is used to map *set* instead. Otherwise this method delegates to :meth:`Homomorphism.image_of_set <thompson.homomorphism.Homomorphism.image_of_set>`.
		
		Examples of finding inverse images:
		
			>>> basis = Generators.standard_basis((2,1))
			>>> basis.expand_to_size(8);
			>>> print(example_5_3.image_of_set(basis, inverse=True))
			[x1 a1 a1 a1 a1, x1 a1 a1 a1 a2 x1 a1 a1 a2 L, x1 a1 a2 a2, x1 a1 a2 a1, x1 a2 a1, x1 a2 a2 a1, x1 a2 a2 a2 a1, x1 a2 a2 a2 a2]
		"""
		if inverse:
			return super().image_of_set(set, self.range.signature, self.domain.signature, self._inv)
		else:
			return super().image_of_set(set)
	
	def repeated_image(self, key, power):
		r"""If :math:`\psi` is the current automorphism, returns :math:`\text{key}\psi^\text{ power}`.
		
		:rtype: a :class:`~thompson.word.Word` instance.
		
		.. doctest::
			
			>>> print(example_4_25.repeated_image('x1 a1', 10))
			x1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1
			>>> print(example_4_25.repeated_image('x1 a1 a1 a1 a1 a1 a1 a1', -3))
			x1 a1
			>>> print(arity_four.repeated_image('x1 a4 a4 a2', 4))
			x1 a3 a3 a2
			>>> print(arity_four.repeated_image('x1 a3 a3 a2', -4))
			x1 a4 a4 a2
		
		.. todo::
			This could be made more efficient for large *power*s by using knowledge of the component containing *key*.
		"""
		inverse = power < 0
		if power == 0:
			if not isinstance(key, Word):
				key = Word(key, self.signature)
			return key
		power = abs(power)
		image = key
		for _ in range(power):
			image = self.image(image, inverse)
		return image
	
	#Group operations
	def __pow__(self, power):
		"""
		>>> psi = random_automorphism()
		>>> ~psi ** 2 == ~psi * ~psi == psi ** -2
		True
		"""
		if not isinstance(power, int):
			raise TypeError("Power must be an integer (instead of {}).".format(power))
		
		if power == 0:
			return Automorphism.identity(self.signature)
		if power == 1:
			return self
		if power == -1:
			return ~self
		
		domain = copy(self.domain)
		range = Generators(self.signature)
		for word in domain:
			range.append(self.repeated_image(word, abs(power)))
		if power > 0:
			pow = Automorphism(domain, range)
			pow.domain_relabeller = self.domain_relabeller
			pow.range_relabeller = self.range_relabeller
		else:
			pow = Automorphism(range, domain)
			pow.domain_relabeller = self.range_relabeller
			pow.range_relabeller  = self.domain_relabeller
		return pow
	
	def __invert__(self):
		"""We overload python's unary negation operator ``~`` as shorthand for inversion. (In Python, ``~`` is normally used for bitwise negation.) We can also call a method explicitily: ``phi.inverse()`` is exactly the same as ``~phi``.
		
			>>> phi = random_automorphism()
			>>> phi * ~phi == ~phi * phi
			True
			>>> (phi * ~phi).is_identity()
			True
			>>> (~phi).quasinormal_basis == phi.quasinormal_basis
			True
		"""
		inv = copy(self)
		inv.domain, inv.range = Generators.sort_mapping_pair(self.range, self.domain)
		inv._map, inv._inv = self._inv, self._map
		
		inv.pond_banks = {(r, -k, ell) for ell, k, r in self.pond_banks}
		inv.characteristics = {(-m, Gamma) for m, Gamma in self.characteristics}
		
		inv.domain_relabeller = self.range_relabeller
		inv.range_relabeller  = self.domain_relabeller
		return inv
	inverse = __invert__
	
	#Quasinormal basis
	#Some ideas for making this process faster:
	#- Cache the component types (wrt QNB) for the QNB elements
	#- Put some part of the core of type B orbits into confirmed? See /theory/Characteristic components and the QNB.tex
	#- Cache the powers (-ell, m) which determine the core part of a component
	def compute_quasinormal_basis(self):
		r"""We say that :math:`\phi` is *in semi-normal form* with respect to the basis :math:`X` if no element of :math:`X` lies in an incomplete :math:`X`-component of a :math:`\phi` orbit. See the :mod:`~thompson.orbits` module for more details.
		
		There is a minimal such basis, :math:`X_\phi` say, and we say that :math:`\phi` is *in quasi-normal form* with respect to :math:`X_\phi`. This method determines and returns the basis :math:`X_\phi` where :math:`\phi` denotes the current automorphism. The result is cached so that further calls to this method perform no additional computation. 
		
			>>> example_4_5.quasinormal_basis
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> alphabet_size_two.quasinormal_basis
			Generators((3, 2), ['x1 a1', 'x1 a2', 'x1 a3', 'x2'])
			>>> example_5_12_phi.quasinormal_basis
			Generators((2, 1), ['x1 a1', 'x1 a2'])
			>>> example_6_2.quasinormal_basis
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2'])
			>>> example_6_8_phi.quasinormal_basis
			Generators((2, 1), ['x1 a1', 'x1 a2'])
		
		:rtype: a :class:`~thompson.generators.Generators` instance.
		
		.. note:: This method additionally looks for ponds and caches the data describing any ponds it finds.
		
		.. todo:: Make the QNB read-only somehow, so that it cannot be expanded once computed.
		
		.. seealso:: Quasi-normal forms are introduced in section 4.2 of the paper. In particular, this method implements Lemma 4.24.1. Higman first described the idea of quasi-normal forms in section 9 of [Hig]_.
		"""
		self.quasinormal_basis = None
		self.pond_banks = None
		
		self.multiplicity = {}
		
		#1. Expand the starting basis until each no element's belongs to a finite X-component.
		basis = self._seminormal_form_start_point()
		basis.cache = set(basis)
		confirmed = set()
		i = 0
		checks_needed = len(basis)
		
		while checks_needed > 0:
			if basis[i] in confirmed:
				i = (i + 1) % len(basis)
				checks_needed -= 1
				continue
			
			ctype, images, _ = self.orbit_type(basis[i], basis)
			if ctype.is_incomplete():
				basis._expand_with_cache(i)
				checks_needed = len(basis)
				continue
			
			elif ctype.is_type_A():
				confirmed.update(images.values())
				#expand basis until each of the images is below basis
				for img in images.values():
					#TODO. I think this operation takes O(len(basis)).
					index, tail = basis.test_above(img, return_index=True)
					for j in range(len(tail)):
						basis._expand_with_cache(index)
						index += -tail[j] - 1
						checks_needed += self.signature.arity - 1
					assert img in basis
				
				#record this new periodic orbit
				period = ctype.characteristic[0]
				try:
					self.multiplicity[period] += 1
				except KeyError:
					self.multiplicity[period] = 1
			
			elif ctype.is_type_B():
				confirmed.add(basis[i])
			i = (i + 1) % len(basis)
			checks_needed -= 1
		
		#Tidy up data
		self.quasinormal_basis = basis
		self.cycle_type = self.multiplicity.keys()
		self.order = float('inf') if len(self.cycle_type) == 0 else lcm(self.cycle_type)
		
		self.characteristics = set()
		for endpt in chain(*self.semi_infinite_end_points()):
			ctype, _, _ = self.orbit_type(endpt, basis)
			if ctype.is_type_B():
				self.characteristics.add(ctype.characteristic)
		
		#2. Look for ponds.
		self.pond_banks = self._find_ponds()
		
		#3. This is naughty. Cast this class as a pure/mixed/infinite aut as appropriate.
		if len(self.multiplicity) == 0:
			from .infinite import InfiniteAut
			self.__class__ = InfiniteAut
		elif len(self.characteristics) == 0:
			from .periodic import PeriodicAut
			self.__class__ = PeriodicAut
		else:
			from .mixed import MixedAut
			self.__class__ = MixedAut
	
	def _find_ponds(self):
		terminal, initial = self.semi_infinite_end_points(exclude_characteristics=True)
		terminal_data = set()
		for term in terminal:
			tail = self._descend_to_complete_infinite(term)
			terminal_data.add((term, tail))
		initial = set(initial)
		
		ponds = []
		for (term, tail) in terminal_data:
			for init in initial:
				power = self._are_banks(term, init, tail)
				if power is not None:
					ponds.append((term, power, init))
					initial.remove(init)
					break
		return ponds
	
	def _seminormal_form_start_point(self):
		r"""Returns the minimal expansion :math:`X` of :math:`\boldsymbol{x}` such that every element of :math:`X` belongs to either *self.domain* or *self.range*. Put differently, this is the minimal expansion of :math:`\boldsymbol{x}` which does not contain any elements which are above :math:`Y \cup W`. See remark 4.10 and example 4.25. This basis that this method produces is the smallest possible which *might* be semi-normal.
		
			>>> example_4_5._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> example_4_11._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1', 'x1 a2'])
			>>> example_4_12._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1', 'x1 a2'])
			>>> example_4_25._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
			>>> cyclic_order_six._seminormal_form_start_point()
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2 a1', 'x1 a1 a2 a2', 'x1 a2'])
		"""
		basis = Generators.standard_basis(self.signature)
		i = 0
		while i < len(basis):
			b = basis[i]
			if b in self.domain or b in self.range:
				i += 1
			else:
				basis.expand(i)
		return basis
	
	def orbit_type(self, y, basis):
		"""Returns the orbit type of *y* with respect to the given *basis*. Also returns a dictionary of computed images, the list (7) from the paper.
		
			>>> #Example 4.5.
			>>> print_component_types(example_4_5, example_4_5.domain)
			x1 a1 a1 a1: Left semi-infinite component with characteristic (-1, a1)
			x1 a1 a1 a2: Bi-infinite component
			x1 a1 a2: Right semi-infinite component with characteristic (1, a2)
			x1 a2 a1: Periodic component of order 2
			x1 a2 a2: Periodic component of order 2
			with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> print_component_types(example_4_5, basis=example_4_5.domain, words=['x', 'x a1', 'x a2'])
			x1: Incomplete component
			x1 a1: Incomplete component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1 a1 a1, x1 a1 a1 a2, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			
			>>> #Example 4.11
			>>> print_component_types(example_4_11)
			x1 a1: Left semi-infinite component with characteristic (-1, a1)
			x1 a2: Right semi-infinite component with characteristic (1, a2)
			with respect to the basis [x1 a1, x1 a2]
		
		.. doctest::
			:hide:
			
			>>> #Example 4.12
			>>> basis = example_4_12._seminormal_form_start_point()
			>>> print_component_types(example_4_12, basis)
			x1 a1: Incomplete component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1, x1 a2]
			>>> basis.expand(0)
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2'])
			>>> print_component_types(example_4_12, basis)
			x1 a1 a1: Bi-infinite component
			x1 a1 a2: Bi-infinite component
			x1 a2: Incomplete component
			with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2]
			>>> basis.expand(2)
			Generators((2, 1), ['x1 a1 a1', 'x1 a1 a2', 'x1 a2 a1', 'x1 a2 a2'])
			>>> print_component_types(example_4_12, basis)
			x1 a1 a1: Periodic component of order 4
			x1 a1 a2: Periodic component of order 4
			x1 a2 a1: Periodic component of order 4
			x1 a2 a2: Periodic component of order 4
			with respect to the basis [x1 a1 a1, x1 a1 a2, x1 a2 a1, x1 a2 a2]
			>>> #Example 4.25
			>>> print_component_types(example_4_25)
			x1 a1: Right semi-infinite component with characteristic (1, a1 a1)
			x1 a2 a1: Bi-infinite component
			x1 a2 a2: Left semi-infinite component with characteristic (-1, a1 a1)
			with respect to the basis [x1 a1, x1 a2 a1, x1 a2 a2]
		
		.. seealso:: Lemmas 4.14, 4.24
		"""
		images = {}
		type_b_data = None
		
		rinf, rpow1, rpow2, rimages = self.test_semi_infinite(y, basis, backward=False)
		for i, image in enumerate(rimages):
			images[i] = image
		
		#Eliminate the periodic case first.
		if rinf and rimages[rpow2] == y:
			assert rpow1 == 0
			ctype = ComponentType.periodic(rpow2)
			del images[rpow2]
			return ctype, images, type_b_data
		
		#Else, we don't have anything periodic.
		linf, lpow1, lpow2, limages = self.test_semi_infinite(y, basis, backward=True)
		for i, image in enumerate(limages):
			images[-i] = image
		
		if not (linf or rinf):
			ctype = ComponentType.incomplete()
			
		elif linf and not rinf:
			tail = y.test_above(limages[-1])
			if tail is None:
				characteristic = None
				type_b, tail = basis.test_above(limages[lpow1])
				type_b_data = type_b_triple(-lpow1, type_b, tail)
			else:
				assert len(tail) > 0
				assert lpow1 == 0
				characteristic = -lpow2, tail
			ctype = ComponentType.semi_infinite(characteristic, backward=True)
		
		elif rinf and not linf:
			tail = y.test_above(rimages[-1])
			if tail is None:
				characteristic = None
				type_b, tail = basis.test_above(rimages[rpow1])
				type_b_data = type_b_triple(rpow1, type_b, tail)
			else:
				assert len(tail) > 0
				assert rpow1 == 0
				characteristic = rpow2, tail
			ctype = ComponentType.semi_infinite(characteristic, backward=False)
		
		elif linf and rinf:
			ctype = ComponentType.complete_infinite()
			type_b, tail = basis.test_above(rimages[rpow1])
			type_b_data = type_b_triple(rpow1, type_b, tail)
		
		#Test_semi_infinite returns one more word than is needed for the 'core'.
		#When appropriate, we remove this extra word.
		#TODO this breaks the conjugacy test and needs to be implemented properly.
		# for index in (rpow2, -lpow2):
			# if abs(index) > 0:
				# try:
					# del images[index]
				# except KeyError:
					# pass
		return ctype, images, type_b_data
	
	def test_semi_infinite(self, y, basis, backward=False):
		r"""Computes the orbit type of *y* under the current automorphism :math:`\psi` with respect to *basis* in the given direction. Let :math:`y\psi^m` be the most recently computed image. The process stops when either:
		
			1. :math:`y\psi^m` is not below the *basis*, for some :math:`m\geq 0`.
			
				- infinite: ``False``
				- start: ``0``
				- end: ``m-1``
				- images: :math:`y, y\psi, \dotsc, y\psi^{m-1}`.
			
			2. For some :math:`0 \le l < m`, :math:`y\psi^l` and :math:`y\psi^m` start with the same basis element.
			
				- infinite: ``True``
				- start: ``l``
				- end: ``m``
				- images: :math:`y, y\psi, \dotsc, y\psi^{m}`.
				
			.. note::
				
				The word :math:`y\psi^m` is not strictly in the core part of the orbit of Lemma 4.24. We return this as part of *images* so that we can compute the characteristic multiplier in :meth:`orbit_type`.
		
		:returns: the tuple *(infinite, start, end, images)*.
		
		.. seealso:: Lemma 4.24
		"""
		image = y
		images = [y]
		m = 0
		heads = set()
		
		result = basis.test_above(y)
		if result is None:
			return False, 0, m-1, images 
		
		prefix, _ = result
		while True:
			heads.add(prefix)
			m += 1
			image = self.image(image, inverse=backward) #Compute y\phi^{\pm m}
			result = basis.test_above(image)
			#1. Did we fall out of X<A>?
			if result is None:
				return False, 0, m-1, images
			
			#2. Otherwise, is there evidence to conclude that this is this an infinite orbit?
			prefix, tail = result
			if prefix in heads:
				for ell, previous in enumerate(images):
					if prefix == previous[:len(prefix)]:
						images.append(image)
						return True, ell, m, images
			images.append(image)
	
	def semi_infinite_end_points(self, exclude_characteristics=False):
		r"""Returns the list of terminal :class:`Words <thompson.word.Word>` in left semi-infinite components and the list of initial words in right semi-infinite components. This is all computed with respect the current automorphism's :meth:`quasinormal_basis`. These are the sets :math:`X\langle A\rangle \setminus Y\langle A\rangle` and :math:`X\langle A\rangle \setminus Z\langle A\rangle`.
		
		:param exclude_characteristics: if True, only the endpoints of non-characteristic semi-infinite orbits will be returned.
		
			>>> print(*example_4_5.semi_infinite_end_points())
			[x1 a1 a1] [x1 a1 a2]
			>>> print(*example_4_11.semi_infinite_end_points())
			[x1 a1] [x1 a2]
			>>> print(*example_4_12.semi_infinite_end_points())
			[] []
			>>> print(*example_4_25.semi_infinite_end_points())
			[x1 a2 a2, x1 a2 a2 a1] [x1 a1, x1 a1 a1]
			>>> print(*example_6_2.semi_infinite_end_points())
			[x1 a1 a1] [x1 a2]
			>>> print(*nathan_pond_example.semi_infinite_end_points())
			[x1 a1 a1, x1 a1 a1 a1, x1 a1 a1 a2] [x1 a1 a2, x1 a1 a2 a1, x1 a1 a2 a2]
			>>> print(*nathan_pond_example.semi_infinite_end_points(exclude_characteristics=True))
			[x1 a1 a1 a2] [x1 a1 a2 a2]
		
		:rtype: A pair of :class:`Generators <thompson.generators.Generators>`.
		
		.. seealso:: The discussion before lemma 4.6.
		"""
		basis = self.quasinormal_basis                  #X
		min_expansion = basis.minimal_expansion_for(self) #Y
		img_expansion = self.image_of_set(min_expansion)  #Z
		terminal = basis.descendants_above(min_expansion) #X<A> \ Y<A>
		initial  = basis.descendants_above(img_expansion) #X<A> \ Z<A>
		
		if exclude_characteristics:
			terminal = Generators(self.signature, (t for t in terminal if not self.orbit_type(t, basis)[0].is_type_B()))
			initial =  Generators(self.signature, (i for i in initial  if not self.orbit_type(i, basis)[0].is_type_B()))
		return terminal, initial
	
	
	def _descend_to_complete_infinite(self, endpt):
		"""A really crude use of the lemma described in AJD's email.
		
		**Lemma.** Let :math:`b` be the bank of a pond. There is a :math:`\Gamma \in A^*` such that :math:`b\Gamma` belongs to a doubly-infinite :math:`\phi`-orbit that meets :math:`X`.
		"""
		basis = self.quasinormal_basis
		for gamma in free_monoid_on_alphas(self.signature.arity):
			w = endpt.extend(gamma)
			ctype, _, _ = self.orbit_type(w, basis)
			if ctype.type == ComponentType._complete_infinite:
				return gamma
	
	def _are_banks(self, term, init, tail):
		#todo docstring
		u = term.extend(tail)
		v = init.extend(tail)
		solns = self.share_orbit(u, v)
		assert not solns.is_sequence(), solns
		if solns.is_empty():
			return None
		power = solns.base
		if self.repeated_image(term, power) == init:
			return power
		else:
			return None
	
	def dump_QNB(self):
		for gen in self.quasinormal_basis:
			ctype, _, _ = self.orbit_type(gen, self.quasinormal_basis)
			print(gen, ctype)
	
	def print_characteristics(self):
		""".. doctest::
			
			>>> example_4_1.print_characteristics()
			(-1, a1)
			(1, a2)
			>>> example_4_25.print_characteristics()
			(-1, a1 a1)
			(1, a1 a1)
			>>> example_6_2.print_characteristics()
			(-2, a1)
			(1, a2)
			>>> (example_6_2 * example_6_2).print_characteristics()
			(-1, a1)
			(1, a2 a2)
			>>> example_6_8_phi.print_characteristics()
			(-1, a1 a1 a1)
			(1, a2 a2 a2)
			>>> #Lemma 5.16
			>>> psi, phi = random_conjugate_pair()
			>>> psi.characteristics == phi.characteristics
			True
		
		.. doctest::
			:hide:
			
			>>> #Lemma 6.1
			>>> from random import randint
			>>> psi = random_infinite_automorphism()
			>>> original_chars = psi.characteristics
			>>> power = psi
			>>> a = randint(2, 6)
			>>> for _ in range(a - 1): 
			... 	power *= psi
			>>> chars = set()
			>>> for mult, char in original_chars:
			... 	d = gcd(mult, a)
			... 	q = abs(a // d)
			... 	chars.add((mult//d, char*q))
			>>> chars == power.characteristics
			True
		
		.. seealso:: Defintion 5.14.
		"""
		for power, mult in sorted(self.characteristics):
			print('({}, {})'.format(power, format(mult)))
	
	#Orbit sharing test
	def share_orbit(self, u, v):
		r"""Determines if :math:`u` and :math:`v` are in the same orbit of the current automorphism :math:`\psi`. Specifically, does there exist an integer :math:`m` such that :math:`u\psi^m = v`?
		
		.. doctest::
			:options: -ELLIPSIS
			
			>>> u  = Word('x1 a2 a3 a1 a2', (3, 2))
			>>> v1 = Word('x1 a1 a2 a2 a3 a1', (3, 2))
			>>> v2 = Word('x2 a3 a2', (3, 2))
			>>> print(alphabet_size_two.share_orbit(u, v1))
			{}
			>>> print(alphabet_size_two.share_orbit(u, v2))
			{-2}
			>>> print(alphabet_size_two.share_orbit(u, u))
			{0}
		
		.. doctest::
			:hide:
			:options: -ELLIPSIS
			
			>>> u  = Word('x a2 a2 a1 a1 a2', (2, 1))
			>>> v1 = Word('x a1 a2', (2, 1))
			>>> v2 = Word('x a1 a1 a2', (2, 1))
			>>> print(example_4_25.share_orbit(u, v1))
			{}
			>>> print(example_4_25.share_orbit(u, v2))
			{3}
			>>> u  = Word('x a2 a2 x a1 a2 L x a2 a1 L x a1 a1 a1 a2 L', (2, 1))
			>>> vs = [
			... 	Word('x a2 a2 a1 a1 a1 a1 a1', (2, 1)),
			... 	Word('x a2 a2 a2', (2, 1)),
			... 	Word('x a2 a2 x a1 a2 L', (2, 1)),
			... 	Word('x a1 a1 x a1 a2 x a2 a2 L L', (2, 1))]
			... 
			>>> for v in vs: print(example_4_25.share_orbit(u, v))
			{-4}
			{}
			{-1}
			{}
			>>> u1 = Word('x a2 a3 a1', (4, 1))
			>>> v1 = Word('x a3 a3 a3 a3 a3 a3 a3 a3 a3 a2 a3 a1', (4, 1))
			>>> v2 = Word('x a1 a2 a3 a4', (4, 1))
			>>> u2 = Word('x a3 a4 a1', (4, 1))
			>>> v3 = Word('x a4 a1 a1', (4, 1))
			>>> v4 = Word('x a4 a3 a2 a1', (4, 1))
			>>> print(arity_four.share_orbit(u1, v1))
			{9}
			>>> print(arity_four.share_orbit(u1, v2))
			{}
			>>> print(arity_four.share_orbit(u2, v3))
			{-1}
			>>> print(arity_four.share_orbit(u2, v4))
			{}
			>>> u = Word("x a1 a2 a1 a2 a1", (2, 1))
			>>> v = Word("x a2 a2 a2 a1", (2, 1))
			>>> print(cyclic_order_six.share_orbit(u, v))
			{..., -1, 2, 5, 8, 11, 14, ...}
			>>> u = Word("x a1 a1 x a1 a2 a1 x a1 a2 a2 a1 L L", (2, 1))
			>>> v1 = Word("x a1 a1 x a2 a2 x a2 a1 L L", (2, 1))
			>>> v2 = Word("x a1 a1 x a1 a2 a1 x a1 a1 L L", (2, 1))
			>>> print(cyclic_order_six.share_orbit(u, v1))
			{..., -1, 5, 11, 17, 23, 29, ...}
			>>> print(cyclic_order_six.share_orbit(u, v2))
			{}
			>>> #Two sides of a pond
			>>> u = Word('x a1 a1 a1 a1 a1 a1 a1 a2', (2, 1))
			>>> v = Word('x a2 a2 a1 a1 a2', (2, 1))
			>>> print(first_pond_example_phi.share_orbit(u, v))
			{4}
		
		:returns: The (possibly empty) :class:`~thompson.orbits.SolutionSet` of all integers :math:`m` for which :math:`u\psi^m = v`. Note that if :math:`u = v` this method returns :math:`\mathbb{Z}`. 
		
		.. seealso:: The implementation is due to lemma 4.24.2 of the paper.
		"""
		#TODO a script which randomly checks examples to verify.
		#TODO pond examples
		orig_u, orig_v = u, v
		basis = self.quasinormal_basis
		if not (basis.is_above(u) and basis.is_above(v)):
			depth = max(u.max_depth_to(basis), v.max_depth_to(basis))
			alphas = range(-1, -self.signature.arity - 1, -1)
			solution_set = SolutionSet.the_integers()
			# For all strings \Gamma of length *depth*:
			for tail in product(alphas, repeat=depth):
				u_desc = u.extend(tail)
				v_desc = v.extend(tail)
				solution_set &= self.share_orbit(u_desc, v_desc)
				if solution_set.is_empty():
					return solution_set
			return solution_set
			#A recursive alternative: intersect and return the results of tests on only the children of u, v.
		
		#Now we're dealing with simple words below the basis.
		u_head, u_tail = basis.test_above(u)
		v_head, v_tail = basis.test_above(v)
		u_head_type, _, u_head_data = self.orbit_type(u_head, basis)
		v_head_type, _, v_head_data = self.orbit_type(v_head, basis)
		
		#Is either element periodic?
		if u_head_type.is_type_A() or v_head_type.is_type_A():
			#If so, are they both peroidic with the same period?
			if u_head_type != v_head_type:
				return SolutionSet.empty_set()
			
			period = u_head_type.characteristic[0]
			image = u
			for i in range(1, period+1):
				image = self.image(image)
				if image == v:
					return SolutionSet(i, period)
			return SolutionSet.empty_set()
		
		#Neither u nor v are periodic.
		#Move along the orbit until we find a nicer u and v.
		u_shift, u_head, u_tail = self._type_b_descendant(u_head, u_tail, u_head_data)
		v_shift, v_head, v_tail = self._type_b_descendant(v_head, v_tail, v_head_data)
		
		u = u_head.extend(u_tail)
		v = v_head.extend(v_tail)
		type, images, _ = self.orbit_type(u, basis)
		
		if self.pond_banks is not None:
			#Add on the other side of a pond if neccessary
			self._check_for_ponds(images)
		for i, image in images.items():
			if image == v:
				assert self.repeated_image(orig_u, u_shift + i - v_shift) == orig_v, (format(orig_u), u_shift + i - v_shift)
				return SolutionSet.singleton(u_shift + i - v_shift)
		return SolutionSet.empty_set()
	
	def _type_b_descendant(self, head, tail, head_data):
		r"""Takes a pair :math:`(y, \Gamma)` below the quasi-normal basis :math:`X` and returns a triple :math:`(\widetilde{y}, \widetilde{\Gamma}, k) where
		* The orbit of :math:`\widetilde{y}` is of type B;
		* :math:`\widetilde{y}\widetilde{\Gamma}` is in the same orbit as :math:`y\Gamma`.
		* :math:`\widetilde{\Gamma}` does not start with the characteristic multiplier of :math:`\widetilde{y}`.
		
			>>> basis = example_4_25.quasinormal_basis
			>>> basis
			Generators((2, 1), ['x1 a1', 'x1 a2 a1', 'x1 a2 a2'])
			>>> head = Word('x a2 a2', (2, 1))
			>>> _, _, type_b_data = example_4_25.orbit_type(head, basis)
			>>> example_4_25._type_b_descendant(head, from_string('a1 a1 a2'), type_b_data)
			(1, Word('x1 a2 a2', (2, 1)), (-2,))
		"""
		orig_head = head
		orig = head.extend(tail)
		if head_data is not None:
			shift_1, head, prefix = head_data['power'], head_data['target'], head_data['end_tail']
			assert self.repeated_image(orig_head, shift_1) == head.extend(prefix), (orig_head, shift_1, head, prefix)
			tail = prefix + tail
		else:
			shift_1 = 0
		
		head_type, _, _ = self.orbit_type(head, self.quasinormal_basis)
		assert head_type.is_type_B(), (head, head_type)
		shift_2, tail = self._remove_from_start(tail, head_type.characteristic)
		assert self.repeated_image(orig, shift_1 + shift_2) == head.extend(tail), (shift_1, shift_2)
		return (shift_1 + shift_2), head, tail
	
	@staticmethod
	def _remove_from_start(tail, characteristic):
		"""Makes a copy of *tail* and removes as many copies of *multiplier* from the start of the copy as possible. Returns a tuple *(k, tail copy)* where *k* is the number of steps forward into the orbit we made.
			
			>>> tail = from_string('a1 a2 a1 a2 a1 a2 a3')
			>>> multiplier = from_string('a1 a2')
			>>> Automorphism._remove_from_start(tail, (2, multiplier))
			(-6, (-3,))
			>>> multiplier = from_string('a1')
			>>> Automorphism._remove_from_start(tail, (2, multiplier))
			(-2, (-2, -1, -2, -1, -2, -3))
		"""
		power, multiplier = characteristic
		i = 0
		n = len(multiplier)
		while tail[i*n : i*n + n] == multiplier:
			i += 1
		return -power*i, tail[i*n:]
	
	def _check_for_ponds(self, images):
		#todo doctests
		M = max(images.keys())
		m = min(images.keys())
		basis = self.quasinormal_basis
		
		for (ell, k, r) in self.pond_banks:
			if images[M] == ell:
				core = self.test_semi_infinite(r, basis)[3]
				for i, img in enumerate(core):
					images[M + k + i] = img
				return
			
			elif images[m] == r:
				core = self.test_semi_infinite(ell, basis, backward=True)[3]
				for i, img in enumerate(core):
					images[m - k - i] = img
				return
	
	def _test_power_conjugate_upto(self, other, sbound, obound, inverses=False):
		r"""In both the periodic and infinite cases, we establish bounds on the powers :math:`a, b` for conjugacy; the rest is brute force.  This method tests to see if :math:`\psi^a` is conjugate to :math:`\phi^b` within the supplied bounds. Should it find a conjugator :math:`\rho`, this method yields a triple :math:`(a, b, \rho)`.
		
		Let :math:`\hat a, \hat b` denote *sbound* and *obound* respectively. If *inverses* is False, then we search over the ranges :math:`1 \le a \le \hat a` and :math:`1 \le b \le \hat b`. If *inverses* is True, we search over the (four times larger) range `1 \le |a| \le \hat a` and :math:`1 \leq |b| \le \hat b`
		"""
		if sbound < 0:
			raise ValueError('sbound parameter should be at least 0 (received {}).'.format(sbound))
		if obound < 0:
			raise ValueError('obound parameter should be at least 0 (received {}).'.format(obound))
		# print(self.__class__.__name__, other.__class__.__name__, sbound, obound)
		if sbound == 0 or obound == 0:
			raise StopIteration
		
		#WLOG let sbound <= obound
		if sbound > obound:
			self, other = other, self
			sbound, obound = obound, sbound
			# print('swapping bounds to', sbound, obound)
			swapped = True
		else:
			swapped = False
		
		s_powers = {1: self}
		if inverses:
			s_powers[-1] = ~self
			inherit_relabellers(s_powers[-1], self)
		for a in range(2, sbound + 1):
			s_powers[a] = s_powers[a-1] * self
			inherit_relabellers(s_powers[a], self)
			if inverses:
				s_powers[-a] = s_powers[1-a] * s_powers[-1]
				inherit_relabellers(s_powers[-a], self)
		
		for b, opow in powers_of(other, obound, inverses):
			for a, spow in s_powers.items():
				# print('trying', a, b)
				# assert spow.domain_relabeller is not None and opow.domain_relabeller is not None
				rho = spow.test_conjugate_to(opow)
				if rho is not None:
					yield (a, b, rho) if not swapped else (b, a, ~rho)

def inherit_relabellers(target, source):
	target.domain_relabeller = source.domain_relabeller
	target.range_relabeller = source.range_relabeller

def type_b_triple(power, head, tail):
	return dict(start_tail = tuple(), power=power, end_tail=tail, target=head)

def powers_of(aut, max, inverses=False):
	assert max >= 1
	yield 1, aut
	power = aut
	for i in range(2, max + 1):
		power *= aut
		inherit_relabellers(power, aut)
		yield i, power
	if not inverses:
		return
	inv = ~aut
	inherit_relabellers(inv, aut)
	yield -1, inv
	power = inv
	for i in range(2, max + 1):
		power *= inv
		inherit_relabellers(power, aut)
		yield -i, power