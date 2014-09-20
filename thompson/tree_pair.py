"""
.. testsetup:: 
	
	from thompson.tree_pair import *
"""
from copy import deepcopy

import svgwrite

from .drawing import *
from .constants import *
from .permutation import Permutation
from .trees import DrawableTree

__all__ = ['TreePair']

class TreePair:
	r"""Thompson's group :math:`V` may be described as a particular set of bijections of the interval :math:`I = [0, 1]`. In turn, these bijections can be represented as pairs of strict binary trees with some additional information. This class is responsible for keeping the three pieces together.
	
	We can form a partition of :math:`I` by slicing it into two halves :math:`[0, 1/2]` and :math:`[1/2, 1]`. If we wanted to, we could slice each of these halves into halves. Equally we could slice some of the newly-formed quarters into eighths, and so on. The result would be a partition
	
	.. math:: 0 = x_0 < x_1 < \dots < x_N < x_{N+1} = 1
	
	where each subinterval :math:`[x_i, x_{i+1}]` looks like :math:`[{a}/{2^n}, (a+1)/{2^n}]` for integers :math:`a` and :math:`n`. Such subintervals are called dyadic; the partition itself is also called dyadic. 
	
	The elements of :math:`V` are the bijections which take two dyadic partitions of :math:`I` and linearly map subintervals from one partition onto subintervals from another. We can use (strict) binary trees to represent dyadic intervals (see :meth:`~thompson.trees.BinaryTree.to_partition`).
	
	The missing ingredient is how the intervals map to each other: for this, we use a :class:`~thompson.permutation.Permutation`.
	
	:ivar int num_leaves: The number of leaves on both trees.
	:ivar tree domain: The domain tree.
	:ivar tree range: The range tree.
	:ivar perm: The :class:`Permutation` of leaves.
		
	"""
	def __init__(self, domain_tree, range_tree, range_labels=None):
		r"""Creates a new tree pair object given a pair of trees and *range_labels*, a string specifying how the leaves of the domain tree map to those of the range tree. Some sanity checks are made on the arguments.
		
		:param domain_tree: a strict DrawableTree, or a string describing one.
		:param range_tree: the same.
		:param range_labels: a string describing the permutation of leaves. If omitted, we assume the leaves' are not permuted at all.
		
		Call the leaves of the domain tree :math:`D_1, D_2, \dotsc, D_N` in depth-first order from left to right---the same order that binary tree :meth:`~thompson.trees.BinaryTree.walk` methods use. For each :math:`i`, label the image in the range tree of :math:`D_i` with an :math:`i`. The argument *range_labels* should be a space-separated string listing the labels of the range tree in depth-first traversal order.
		
		.. figure:: examples/tree_pair/TreePair_render.svg
			
			**Example.** This is the object ``TreePair("11000", "10100", "1 2 3")``.
		
		:raises ValueError:
		- if the two trees have a different number of leaves,
		- if the trees are given (i.e. not generated from a string) and one of them is not :meth:`strictly binary <thompson.trees.BinaryTree.is_strict>`,
		- if range_labels doesn't properly describe a :class:`~thompson.permutation.Permutation`,
		- if range_labels describes a permutation of too many/few leaves.
		"""
		if isinstance(domain_tree, str):
			domain_tree = DrawableTree.from_string(domain_tree)
		elif not domain_tree.is_strict():
			raise ValueError("Domain tree is not strictly binary.")
		
		if isinstance(range_tree, str):
			range_tree = DrawableTree.from_string(range_tree)
		elif not domain_tree.is_strict():
			raise ValueError("Range tree is not strictly binary.")
		
		if domain_tree.num_leaves() != range_tree.num_leaves():
			raise ValueError("Domain tree has %i leaves, but range tree has %i leaves." 
			  % (domain_tree.num_leaves(), range_tree.num_leaves()))
		
		self.num_leaves = domain_tree.num_leaves()
		self.domain = domain_tree
		self.range = range_tree
		
		if range_labels is None:
			range_labels = list(range(1, self.num_leaves+1))
		self.perm = Permutation(range_labels).inverse()
		if self.perm.size != self.num_leaves:
			raise ValueError("range_labels permutes %i leaves, but the trees have %i leaves."
			  % (self.num_leaves, self.perm.size))
	
	@creates_svg
	def render(self, name=None, **kwargs):
		"""Renders a representation of the group element in terms of trees. The domain and range trees are rendered, and an arrow (with optional label *name*) is drawn between them. Labels are added to the leaves to describe the permutation.
		
		.. figure:: examples/tree_pair/TreePair_render.svg
			
			**Example.** The output of ``TreePair("11000", "10100", "1 2 3").render()``. [:download:`Source code <examples/tree_pair/TreePair_render.py>`].
		"""
		#1. Setup. Create a container group and get the SVG groups for the trees.
		container = svgwrite.container.Group(class_='element')
		
		self.domain.layout()
		self.range.layout()
		
		left = self.domain.render(leaf_names=range(1, self.num_leaves + 1))
		right = self.range.render(leaf_names=self.perm.inverse().output)
		
		left['class'] = 'domain'
		right['class'] = 'range'
		
		container.add(left)
		container.add(right)
		
		#2. Draw an arrow between the two trees.
		if name is not None:
			mid = Coord.unscaled(left.size.x, min(left.size.y, right.size.y)/2) + Coord(ARROW_LENGTH/2, 0)
			start = Coord(-ARROW_LENGTH/2, 0)
			end = -start

			arrow_parent = svgwrite.container.Group()
			arrow_parent.translate(mid)
			label = svgwrite.text.Text(name, (0, 0), class_="above")
			arrow_parent.add(label)
			container.add(arrow_parent)
		else:
			start = left.size.scale(1, 0.5)
			end = start + Coord(ARROW_LENGTH, 0)
			arrow_parent = container
		
		arrow = svgwrite.shapes.Line(start, end, class_='arrow')
		arrow_parent.add(arrow)
		
		#3. Finally, position the right tree.
		start = left.size.to_x() + Coord(ARROW_LENGTH, 0)
		right.translate(start)
		
		size = Coord(start.x + right.size.x, max(left.size.y, right.size.y), scale=1)
		set_size(container, size)
		return container
	
	@creates_svg
	def render_bijection(self, **kwargs):
		"""Returns an SVG group containing a plot of *self*, rendered as a bijection of :math:`[0, 1]`.
		
		.. figure:: examples/tree_pair/TreePair_render_bijection.svg
			
			**Example.** The output of ``TreePair("11000", "10100", "1 2 3").render_bijection()``. [:download:`Source code <examples/tree_pair/TreePair_render.py>`].
		"""
		#0. Setup.
		g = svgwrite.container.Group(class_='graph')
		
		x_partition, _ , _ = self.domain.to_partition()
		y_partition, _ , _ = self.range.to_partition()
				
		#1. Draw both the axes.
		x_axis = svgwrite.shapes.Polyline(class_='axis')
		for x in x_partition:
			mark = Coord(x, 0) * GRAPH_SCALE_FACTOR
			x_axis.points.append(mark)
			label = svg_fraction(x, Coord(0, 0.5) + mark)
			g.add(label)
		g.add(x_axis)
		
		y_axis = svgwrite.shapes.Polyline(class_='axis')
		for y in y_partition:
			mark = Coord(0, -y) * GRAPH_SCALE_FACTOR
			y_axis.points.append(mark)
			label = svg_fraction(y, Coord(-0.5, 0) + mark)
			g.add(label)
		g.add(y_axis)
		
		#2. Draw the segments of the plotted function.
		start = end = None
		for i in range(1, self.num_leaves+1):
			#Draw the ith segment of the graph.
			x = x_partition[i-1]
			y = y_partition[self.perm[i]-1]
			start = Coord(x, -y)*GRAPH_SCALE_FACTOR
			
			new_segment = start != end
			x = x_partition[i]
			y = y_partition[self.perm[i]]
			end = Coord(x, -y) * GRAPH_SCALE_FACTOR
			
			if new_segment:
				segment = svgwrite.shapes.Polyline(points=(start,), class_='plot')
				g.add(segment)
			segment.points.append(end)
		
		size = Coord(1, 1) * GRAPH_SCALE_FACTOR + Coord(2, 2)
		set_size(g, size, offset=Coord(0, GRAPH_SCALE_FACTOR) + Coord(1, 1))
		return g
	
	def in_F(self):
		"""Returns True if this represents an element of *F*. This happens when no permutation of the leaves occurs.
		
			>>> TreePair("11000", "10100", "1 2 3").in_F() #identity permutation
			True
			>>> TreePair("11000", "10100", "2 3 1").in_F() #3-cycle (1 2 3)
			False
			>>> TreePair("11000", "10100", "1 3 2").in_F() #2-cycle (2 3)
			False
		"""
		return self.perm.is_identity()
	
	def in_T(self):
		"""Returns True if this represents an element of *T*. This happens when all the leaves are shifted by the same amount (possibly zero).
		
			>>> TreePair("11000", "10100", "1 2 3").in_T() #identity permutation
			True
			>>> TreePair("11000", "10100", "2 3 1").in_T() #3-cycle (1 2 3)
			True
			>>> TreePair("11000", "10100", "1 3 2").in_F() #2-cycle (2 3)
			False
		"""
		shift = self.perm.output.index(1)
		for index, value in self.perm:
			if value % self.perm.size != (index - shift) % self.perm.size:
				return False
		return True
	
	def reduce(self):
		"""Many different tree pairs correspond to the same bijection of :math:`[0, 1]`. We can sometimes remove nodes from a tree pair to yield a pair of smaller trees which represents the same  bijection. We remove matching pairs of carets until it is no longer possible to do so. At this point, we call the tree pair *reduced*. Cannon, Floyd and Parry show in section 2 of [CFP]_ that there is precisely one reduced tree diagram for each element of F.
		
		This method *modifies* the TreePair so that it becomes reduced. The examples below create an unreduced and reduced tree, then check to see if ``unreduced`` == ``reduced``. The :meth:`equality comparison <__eq__>` reduces both trees to perform the comparison.
		
			>>> w = TreePair("111100000", "111100000") #Identity element
			>>> w == TreePair("0","0")
			True
			>>> x = TreePair("1110010010100", "1111100000100") #in F
			>>> x == TreePair("110100100", "111100000")
			True
			>>> y = TreePair("11100100100", "11100010100", "3 4 5 6 1 2") #in T
			>>> y == TreePair("1100100", "1100100", "2 3 4 1")
			True
			>>> z = TreePair("110110010010100", "111001001100100", "2 3 7 8 1 6 4 5") #in V
			>>> z == TreePair("110100100", "110011000", "2 5 1 4 3")
			True
		"""
		d_leaves = self.domain.leaves()
		#TODO Maybe the method that does this should belong to the permutation class?
		r_leaves = self.range.leaves(perm=self.perm)
		i = 0
		#i   = 0, 1, ...., size-2 (pair indices)
		#i+1 = 1, ...., size-1    (permutation indices)
		while i < self.num_leaves - 1:
			#If the (i, i+1)th leaves form a caret on the domain tree,
			#and if their images form a caret on the range tree
			#and if the left/right order is preserved,
			if (d_leaves[i].parent is d_leaves[i+1].parent and
			  r_leaves[i].parent is r_leaves[i+1].parent and
			  r_leaves[i].is_left_child()):
				#then both carets can be removed.
				self._delete_caret(d_leaves, i)
				self._delete_caret(r_leaves, i)
				self.perm.remove_image_of(i+2)
				self.num_leaves -= 1
				#Need to backtrack in case this contraction has formed a new caret behind it.
				if i > 0:
					i -= 1
			else:
				i += 1
	
	@staticmethod
	def _delete_caret(list, index):
		del list[index + 1]  
		list[index] = list[index].parent
		list[index].detach_children()
	
	def __eq__(self, other):
		"""Two TreePairs are equal if they have the same permutation, domain trees and range trees after reduction. See :meth:`reduce()` for examples.
		"""
		
		if not isinstance(other, TreePair):
			return NotImplemented
		self.reduce()
		other.reduce()
		return (self.perm   == other.perm   and
		        self.domain == other.domain and 
		        self.range  == other.range)
		
	def __ne__(self, other):
		if not isinstance(other, TreePair):
			return NotImplemented
		return not (self == other)
	
	def is_identity(self):
		""":meth:`reduce` s this tree pair. Returns true if the pair represents the identity element of V.
		
			>>> TreePair("10100", "10100").is_identity() #identical trees, no permutation
			True
			>>> TreePair("111100000", "111100000").is_identity() #identical trees, no permutation
			True
			>>> TreePair("111100000", "111100000", "2 3 4 5 1").is_identity() #permute leaves
			False
		"""
		self.reduce()
		return self.domain.is_trivial() and self.range.is_trivial()
	
	def __mul__(self, other): #self * other
		r"""We use Python's multiplication to represent composition of tree pairs. The order is backwards compared to the usual meaning of composition. If ``F`` and ``G`` are tree pairs representing functions :math:`f` and :math:`g\colon` :math:`[0, 1] \to [0, 1]` , then ``F * G`` represents the composition
		
			.. math:: g \circ f : t \mapsto g(f(t))
		
		('f, then g' rather than 'f after g').
		
		**NB.** Both *self* and *other* are reduced before the multiplication begins. The product is also reduced before it is returned.
		
			>>> #Trivial * trivial = trivial
			>>> TreePair("0", "0") * TreePair("100", "100") == TreePair("0", "0")
			True
			>>> #Same trees, two different permutations
			>>> x = TreePair("1100100", "1100100", "1 4 2 3") * TreePair("1100100", "1100100", "2 3 4 1")
			>>> x == TreePair("1100100", "1100100", "4 2 3 1")
			True
			>>> #Same trees within each pair
			>>> y = TreePair("100", "100", "2 1") * TreePair("10100", "10100", "2 3 1")
			>>> y == TreePair("11000", "10100", "3 2 1")
			True
			>>> #Complicated trees, identity permutations
			>>> z = TreePair("1101000", "1100100") * TreePair("1110000", "1011000")
			>>> z == TreePair("111001000", "101100100")
			True
		"""
		if not isinstance(other, TreePair):
			return NotImplemented
		
		self.reduce()
		other.reduce()
		
		s = deepcopy(self)
		o = deepcopy(other)
		
		s.expand(o)
		# assert s.range == o.domain, "Trees not equal"
		prod = TreePair(s.domain, o.range)
		# print(s.perm, o.perm)
		prod.perm = s.perm * o.perm
		return prod
	
	def expand(self, other, sran=None, odom=None, s_inserted = 0, sdom_leaves=None, o_inserted = 0, oran_leaves=None):
		"""Expands two tree pairs so that they can be multiplied."""
		if sran is None: sran = self.range
		if odom is None: odom = other.domain
		if sdom_leaves is None: sdom_leaves = self.domain.leaves(perm=self.perm.inverse())
		if oran_leaves is None: oran_leaves = other.range.leaves(perm=other.perm) #probly
		
		# print('expand:', sran.name, odom.name, s_inserted, names(sdom_leaves), o_inserted, names(oran_leaves))
		
		if sran.is_leaf() and not odom.is_leaf():
			#replace the preimage of sran by a copy of odom
			# print('copy', odom.name, 'onto', sdom_leaves[0].name)
			
			subtree = deepcopy(odom)
			# subtree.name = "copy of" + subtree.name
			sdom_leaves[0].replace_with(subtree)
			
			sdom_leaves.pop(0)# print('Removing', name(), 'from sdom_leaves')
			for child in subtree.walk():
				if child.is_leaf(): 
					oran_leaves.pop(0)# print("subtree contains leaf", name(), "removing from oran_leaves")
			
			insertion_count = subtree.num_leaves() - 1
			s_inserted += insertion_count
			self.num_leaves += insertion_count
			# print('inserted', insertion_count, 'leaves to self')
			
			#update the permutation
			# print("expanding self permutation index", s_inserted+1, "to width", insertion_count + 1)
			# print("from", repr(self.perm), end=" ")
			self.perm.expand_range(s_inserted+1, insertion_count + 1)
			# print('to', repr(self.perm))
		
			
		elif not sran.is_leaf() and odom.is_leaf():
			#replace the image of odom by sran
			# print('copy', sran.name, 'onto', oran_leaves[0].name)
			subtree = deepcopy(sran)
			# subtree.name = "copy of" + subtree.name
			insertion_count = subtree.num_leaves() - 1
			oran_leaves[0].replace_with(subtree)
			
			oran_leaves.pop(0)# print('Removing', name(), 'from oran_leaves')
			for child in subtree.walk():
				if child.is_leaf():
					sdom_leaves.pop(0)# print("subtree contains leaf", name(), "removing from sdom_leaves")
			
			insertion_count = subtree.num_leaves() - 1
			o_inserted += insertion_count
			self.num_leaves += insertion_count
			# print('inserted', insertion_count, 'leaves to other')
			
			#update the permutation
			# print("expanding other permutation from", repr(other.perm), end=" ")
			other.perm.expand_domain(o_inserted+1, insertion_count + 1)
			# print('to', repr(other.perm))
			
		elif not sran.is_leaf() and not odom.is_leaf():
			s_inserted, sdom_leaves, o_inserted, oran_leaves =\
			    self.expand(other, sran.left, odom.left, s_inserted, sdom_leaves, o_inserted, oran_leaves)
			s_inserted, sdom_leaves, o_inserted, oran_leaves =\
			    self.expand(other, sran.right, odom.right, s_inserted, sdom_leaves, o_inserted, oran_leaves)
		
		elif sran.is_leaf() and odom.is_leaf():
			sdom_leaves.pop(0)# print('Removing', name(), 'from sdom_leaves')
			oran_leaves.pop(0)# print('Removing', name(), 'from oran_leaves')
		
		# print('returning', s_inserted, o_inserted, names(sdom_leaves), names(oran_leaves))
		return s_inserted, sdom_leaves, o_inserted, oran_leaves
		
#TODO: compose, invert

def name(x):
	if hasattr(x, 'name'):
		return x.name
	return repr(x)

def names(xs):
	return [name(x) for x in xs]

#Named TreePairs
#TODO check these
#TODO label these as constants - not to be modified?
A = TreePair("10100", "11000")
B = TreePair("1010100", "1011000")
C = TreePair("10100", "10100", "2 3 1")
pi_0 = TreePair("10100", "10100", "2 1 3")

#TODO: memoize?
def x(n):
	"""todo check me"""
	return TreePair("10" * n + "100", "10" * (n-1) + "11000")