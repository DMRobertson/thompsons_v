from .trees import DrawableTree
from .drawing import *

__all__ = ['TreePair']

class TreePair:
	r"""Thompson's group :math:`V` may be described as a particular set of bijections of the interval :math:`I = [0, 1]`. In turn, these bijections can be represented as pairs of strict binary trees with some additional information. This class is responsible for keeping the three pieces together.
	
	We can form a partition of :math:`I` by slicing it into two halves :math:`[0, 1/2]` and :math:`[1/2, 1]`. If we wanted to, we could slice each of these halves into halves. Equally we could slice some of the newly-formed quarters into eighths, and so on. The result would be a partition
	
	.. math:: 0 = x_0 < x_1 < \dots < x_N < x_{N+1} = 1
	
	where each subinterval :math:`[x_i, x_{i+1}]` looks like :math:`[{a}/{2^n}, (a+1)/{2^n}]` for integers :math:`a` and :math:`n`. Such subintervals are called dyadic; the partition itself is also called dyadic. 
	
	The elements of :math:`V` are the bijections which take two dyadic partitions of :math:`I` and linearly map subintervals from one partition onto subintervals from another. We can use (strict) binary trees to represent dyadic intervals (see :meth:`~thompson.trees.BinaryTrees.to_partition`). The missing ingredient is how the intervals map to each other.
	"""
	def __init__(self, domain_tree, range_tree, range_labels):
		r"""Creates a new tree pair object given a pair of trees and *range_labels*, a string specifying how the leaves of the domain tree map to those of the range tree. Some sanity checks are made on the arguments.
		
		:param domain_tree: a strict DrawableTree, or a string describing one.
		:param range_tree: the same.
		:param range_labels: a string specifying how the leaves are mapped from the domain to range tree.
		
		Label the leaves of the domain tree :math:`D_1, D_2,\dotsc, D_N` in depth-first order from left to right---the same order that binary tree :meth:`~thompson.trees.BinaryTree.walk` methods use. For each :math:`i`, label the image in the range tree of :math:`D_i` with an :math:`i`. The argument *range_labels* should be a space-separated string listing the labels of the range tree in depth-first traversal order.
		
		.. figure:: examples/example_permutation.png
			
			**Example.** This is the object ``TreePair("1100100", "1110000", "1 2 0 3")``.
		
		:raises ValueError: if the two trees have a different number of leaves.
		
		.. todo::
			
		"""
		if isinstance(domain_tree, str):
			domain_tree = DrawableTree.from_string(domain_tree)
		if isinstance(range_tree, str):
			range_tree = DrawableTree.from_string(range_tree)
		
		if domain_tree.num_leaves() != range_tree.num_leaves():
			raise ValueError("Domain tree has %i leaves, but range tree has %i leaves." \
			  % (self.num_leaves, range_tree.num_leaves()))
		
		self.num_leaves = self.domain.num_leaves()
		self.domain = domain_tree
		self.range = range_tree
		
		self.range_labels = [int(x) for x in range_labels.split()]
		assert len(self.range_labels) == self.num_leaves, "Permutation is not fully specified."
		
		self.perm = {}
		for traversal_index, label in enumerate(self.range_labels):
			self.perm[label] = i
	
	@creates_SVG
	def render(self, name='test'):
		#1. Create the drawing and render the two trees.
		dwg, canvas = new_drawing(filename, True)
		y_offset = self.layout(canvas, name)
		
		"""if plot_bijection:
			graph = self.plot_bijection()
			graph.translate(Coord(1, y_offset + GRAPH_SCALE_FACTOR))
			canvas.add(graph)"""
		
	@creates_SVG
	def layout(self, name):
		"""Renders the two trees and positions them in a group."""
		#1. Setup. Create a container group and get the SVG groups for the trees.
		container = svgwrite.container.Group(class_='element')
		
		self.domain.layout()
		self.range.layout()
		
		left = self.domain.render(leaf_names=range(self.num_leaves))
		right = self.range.render(leaf_names=self.range_labels)
		
		left['class'] = 'domain'
		right['class'] = 'range'
		
		container.add(left)
		container.add(right)
		
		#2. Position the two trees and draw an arrow between them.
		w = left.size[0]
		y = (left.size[1] - 1)/2
		right.translate(Coord(w + ARROW_LENGTH + 1, 0))
		
		g = svgwrite.container.Group(class_='group_element')
		canvas.add(g)
		
		mid = coord(w + ARROW_LENGTH/2, y)
		start = coord(-ARROW_LENGTH/2, 0)
		end   = coord(+ARROW_LENGTH/2, 0)
		g.translate(mid)
		
		arrow = svgwrite.shapes.Line(start, end, class_='arrow')
		g.add(arrow)
		
		if name is not None:
			g.add(svgwrite.text.Text(name, (0, 0), class_="above"))
		
		y_offset = max(left.size[1], right.size[1]) + 1
		return container
	
	@creates_SVG
	def render_bijection(self):
		"""Returns an SVG group."""
		#3. Plot the group element as a set bijection of [0, 1]
		g = svgwrite.container.Group(class_='graph')
		
		x_partition, _ , _ = self.domain.to_partition()
		y_partition, _ , _ = self.range.to_partition()
		
		assert len(x_partition) == len(y_partition) == self.num_leaves + 1, "Partitions lengths improper."
		
		x_axis = svgwrite.shapes.Polyline(class_='axis')
		for x in x_partition:
			mark = Coord(x, 0, scale=GRAPH_SCALE_FACTOR)
			x_axis.points.append(mark)
			g.add(svgwrite.text.Text(x, insert=mark, class_='below'))
		g.add(x_axis)
		
		y_axis = svgwrite.shapes.Polyline(class_='axis')
		for y in y_partition:
			mark = Coord(0, -y, scale=GRAPH_SCALE_FACTOR)
			y_axis.points.append(mark)
			g.add(svgwrite.text.Text(str(y) + " ", insert=(-1*ex, mark[1]), class_='left centered'))
		g.add(y_axis)
		
		start = end = None
		for i in range(self.num_leaves):
			#Draw the ith segment of the graph.
			x = x_partition[i]
			y = y_partition[self.perm[i]]
			start = Coord(x, -y, scale=GRAPH_SCALE_FACTOR)
			
			new_segment = start != end
			
			x = x_partition[i+1]
			y = y_partition[self.perm[i]+1]
			end   = Coord(x, -y, scale=GRAPH_SCALE_FACTOR)
			
			if new_segment:
				segment = svgwrite.path.Path(d=('M', start))
				g.add(segment)
			segment.push('L', end)
		
		return g