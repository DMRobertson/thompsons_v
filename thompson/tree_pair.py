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
		
		:raises ValueError: if the two trees have a different number of leaves,
		:raises ValueError: if the trees are given (i.e. not generated from a string) and one of them is not :meth:`strictly binary <thompson.trees.BinaryTree.is_strict>`
		:raises ValueError: if range_labels doesn't properly describe a :class:`~thompson.permutation.Permutation`,
		:raises ValueError: if range_labels describes a permutation of too many/few leaves.
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
		
		x_partition = self.domain.to_partition()
		y_partition = self.range.to_partition()
				
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
		
		**NB.** Both *self* and *other* are reduced before the multiplication begins. TODO: The product is also reduced before it is returned.
		
			>>> #Trivial * trivial = trivial
			>>> TreePair("0", "0") * TreePair("100", "100") == TreePair("0", "0")
			True
			>>> #Same trees, two different permutations
			>>> x = TreePair("1100100", "1100100", "1 4 2 3") * TreePair("1100100", "1100100", "2 3 4 1")
			>>> x == TreePair("1100100", "1100100", "4 2 3 1")
			True
			>>> #Same trees within each pair
			>>> y = TreePair("100", "100", "2 1") * TreePair("10100", "10100", "3 1 2")
			>>> y == TreePair("11000", "10100", "2 3 1")
			True
			>>> #Complicated trees, identity permutations
			>>> z = TreePair("1101000", "1100100") * TreePair("1110000", "1011000")
			>>> z == TreePair("111001000", "101100100")
			True
			>>> #Powers of a permutation
			>>> a = TreePair("11010100100", "11010100100", "1 6 5 2 3 4")
			>>> a * a * a * a * a * a == TreePair("0", "0")
			True
		
		.. _multiplication example: :download:
		
		.. figure:: examples/tree_pair/TreePair_mul.svg
			:alt: Drawings of f, g, and the compositions fg and gf.
			:target: _downloads/TreePair_mul.svg
		
			**Example.** Let *f* be ``TreePair("100", "100", "2 1")`` and let *g* be ``TreePair("10100", "10100", "1 3 2")``. In terms of tree pairs, *f* exchanges the subtrees of the root; *g* does the same to the right child of the root. There are two ways to compose these elements, yielding two (in this case different) products.
			
			[:download:`Source code <examples/tree_pair/TreePair_mul.py>`] [:download:`Full size <examples/tree_pair/TreePair_mul.svg>`].
		
		
		"""
		# TODO: really horrible example (4 different trees and 4 messy permutations)
		# TODO: TreePair.__str__() and __repr__()
		if not isinstance(other, TreePair):
			return NotImplemented
		
		self.reduce()
		other.reduce()
		
		s = deepcopy(self)
		o = deepcopy(other)
		s._expand(o)
		
		prod = TreePair(s.domain, o.range)
		prod.perm = o.perm * s.perm #o after s
		#TODO. Is prod reduced at this stage? Should it be? Hmm.
		return prod
	
	def _expand(self, other, sran=None, odom=None, s_inserted = 0, sdom_leaves=None, o_inserted = 0, oran_leaves=None):
		"""Expands two tree pairs so that they can be multiplied. The general idea:
		
		1. Let sran = self.range, and odom=other.domain.
		2. Are both sran and odom branches? If so, call this function again with:
			a. sran = sran.left, odom = odom.left
			b. sran = sran.right, odom = odom.right
		3. Else, if sran is a leaf and odom is not:
			a. Copy odom onto the preimage under self of sran.
			b. Update the permutation of self.
		4. Else, if odom is a leaf and sran is not:
			a. Copy sran onto the image under other of sdom.
			b. Update the permutation of other.
		5. Else, both nodes are branches. Do nothing.
		6. Return variables describing the current traversal state.
		"""
		if sran is None: sran = self.range;
		if odom is None: odom = other.domain;
		if sdom_leaves is None: sdom_leaves = self.domain.leaves(perm=self.perm.inverse())
		if oran_leaves is None: oran_leaves = other.range.leaves(perm=other.perm) #probly
		
		#A module constant DEBUG_MULTIPLICATION toggles the various debug messages
		if not DEBUG_MULTIPLICATION:
			global print
			print = no_op
		else:
			#Assign coordinates for print-out messages
			self.domain.layout()
			self.range.layout()
			other.domain.layout()
			other.range.layout()
		print('expanding sran ({},{}) and odom ({},{}).'.format(sran.x, sran.y, odom.x, odom.y),
			"Different:", sran == odom)
		print("sdom_leaves", [(l.x, l.y) for l in sdom_leaves])
		
		if sran.is_leaf() and not odom.is_leaf():
			# replace the preimage of sran by a copy of odom
			print('copy odom onto sdom')
			subtree = deepcopy(odom)
			sdom_leaves[0].replace_with(subtree)
			sdom_leaves.pop(0)
			print("sdom_leaves pop")
			
			#Mark any leaves we copy over as dealt with
			for i in range(subtree.num_leaves()):
				oran_leaves.pop(0)
				print("oran_leaves pop")
			
			insertion_count = subtree.num_leaves() - 1
			s_inserted += insertion_count
			self.num_leaves += insertion_count
			print('inserted', insertion_count, 'leaves to self')
			
			#update the permutation
			index = self.num_leaves - s_inserted - len(sdom_leaves)
			print("expanding self permutation range", index, "to width", insertion_count + 1)
			print("from", repr(self.perm), end=" ")
			self.perm.expand_range(index, insertion_count + 1)
			print('to', repr(self.perm))
		
		elif (not sran.is_leaf()) and odom.is_leaf():
			#replace the image of odom by sran
			print('copying sran onto oran')
			subtree = deepcopy(sran)
			oran_leaves[0].replace_with(subtree)
			oran_leaves.pop(0)
			print("oran_leaves pop")
			
			for i in range(subtree.num_leaves()):
				sdom_leaves.pop(0)
				print("sdom_leaves pop")
			
			insertion_count = subtree.num_leaves() - 1
			o_inserted += insertion_count
			other.num_leaves += insertion_count
			print('inserted', insertion_count, 'leaves to other')
			
			#update the permutation
			index = other.num_leaves - o_inserted - len(oran_leaves)
			print("expanding other permutation index", index, "to width", insertion_count + 1)
			print("from", repr(other.perm), end=" ")
			other.perm.expand_domain(index, insertion_count + 1)
			print('to', repr(other.perm))
			
		elif (not sran.is_leaf()) and not odom.is_leaf():
			# print('left child')
			s_inserted, sdom_leaves, o_inserted, oran_leaves =\
			    self._expand(other, sran.left, odom.left, s_inserted, sdom_leaves, o_inserted, oran_leaves)
			# print('right child')
			s_inserted, sdom_leaves, o_inserted, oran_leaves =\
			    self._expand(other, sran.right, odom.right, s_inserted, sdom_leaves, o_inserted, oran_leaves)
		
		elif sran.is_leaf() and odom.is_leaf():
			print('both sran and odom are leaves leaves')
			sdom_leaves.pop(0)
			oran_leaves.pop(0)
			print("sdom_leaves pop")
			print("oran_leaves pop")
		print('returning to parent')
		return s_inserted, sdom_leaves, o_inserted, oran_leaves
	
	def inverse(self):
		"""Copies this tree pair, then modifies the copy so that it's the inverse of the original pair. The inverse is returned.
		
			>>> TreePair("11000", "10100", "2 3 1").inverse()
			TreePair('10100', '11000', '3 1 2')
		"""
		x = deepcopy(self)
		x.domain, x.range = x.range, x.domain
		x.perm.invert()
		return x

DEBUG_MULTIPLICATION = False
"""Set this to ``True`` to print debug messages when :meth:`multiplying <TreePair.__mul__>` two tree pairs."""

def no_op(*args, **kwargs):
	pass

#TODO: random tree pair

#Named TreePairs
#TODO check these
#TODO label these as constants - not to be modified?
A = TreePair("10100", "11000")
B = TreePair("1010100", "1011000")
C = TreePair("10100", "10100", "2 3 1")
pi_0 = TreePair("10100", "10100", "2 1 3")

# TODO loads more examples

#TODO: memoize?
def x(n):
	"""todo check me"""
	return TreePair("10" * n + "100", "10" * (n-1) + "11000")