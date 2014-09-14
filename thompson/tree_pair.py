"""
.. testsetup:: 
	
	from thompson.tree_pair import *
"""

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
	
	The elements of :math:`V` are the bijections which take two dyadic partitions of :math:`I` and linearly map subintervals from one partition onto subintervals from another. We can use (strict) binary trees to represent dyadic intervals (see :meth:`~thompson.trees.BinaryTree.to_partition`). The missing ingredient is how the intervals map to each other.
	
	:ivar int num_leaves: The number of leaves on both trees.
	:ivar tree domain: The domain tree.
	:ivar tree range: The range tree.
	:ivar perm: The :class:`Permutation` of leaves.
		
	"""
	def __init__(self, domain_tree, range_tree, range_labels):
		r"""Creates a new tree pair object given a pair of trees and *range_labels*, a string specifying how the leaves of the domain tree map to those of the range tree. Some sanity checks are made on the arguments.
		
		:param domain_tree: a strict DrawableTree, or a string describing one.
		:param range_tree: the same.
		:param range_labels: a string specifying how the leaves are mapped from the domain to range tree.
		
		Call the leaves of the domain tree :math:`D_1, D_2, \dotsc, D_N` in depth-first order from left to right---the same order that binary tree :meth:`~thompson.trees.BinaryTree.walk` methods use. For each :math:`i`, label the image in the range tree of :math:`D_i` with an :math:`i`. The argument *range_labels* should be a space-separated string listing the labels of the range tree in depth-first traversal order.
		
		.. figure:: examples/tree_pair/TreePair_render.svg
			
			**Example.** This is the object ``TreePair("11000", "10100", "1 2 3")``.
		
		:raises ValueError: if the two trees have a different number of leaves; if range_labels doesn't properly describe a permutation (see :class:`~thompson.permutation.Permutation`); or if range_labels describes a permutation of too many/few leaves.
		"""
		if isinstance(domain_tree, str):
			domain_tree = DrawableTree.from_string(domain_tree)
		if isinstance(range_tree, str):
			range_tree = DrawableTree.from_string(range_tree)
		
		if domain_tree.num_leaves() != range_tree.num_leaves():
			raise ValueError("Domain tree has %i leaves, but range tree has %i leaves." 
			  % (self.num_leaves, range_tree.num_leaves()))
		
		self.num_leaves = domain_tree.num_leaves()
		self.domain = domain_tree
		self.range = range_tree
		
		self.perm = Permutation(range_labels)
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
		right = self.range.render(leaf_names=self.perm.output)
		
		left['class'] = 'domain'
		right['class'] = 'range'
		
		container.add(left)
		container.add(right)
		
		#2. Draw an arrow between the two trees.
		if name is not None:
			mid = left.size.scale(1, 0.5) + Coord(ARROW_LENGTH/2, 0)
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
			
			**Example.** The output of ``TreePair("11000", "10100", "1 2 3").render_bijection)``. [:download:`Source code <examples/tree_pair/TreePair_render.py>`].
		"""
		#0. Setup.
		g = svgwrite.container.Group(class_='graph')
		
		x_partition, _ , _ = self.domain.to_partition()
		y_partition, _ , _ = self.range.to_partition()
		
		#TODO make the trees perform this check
		assert len(x_partition) == len(y_partition) == self.num_leaves + 1, "Partitions lengths improper."
		
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
		return self.perm.is_identity() or self.perm.is_cycle(of_length=self.num_leaves)
