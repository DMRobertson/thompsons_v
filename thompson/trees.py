"""In graph theory, the word `tree <http://en.wikipedia.org/wiki/Tree_(graph_theory)>`__ means a connected graph with no cycles. In computer science (CS), a `tree <http://en.wikipedia.org/wiki/Tree_structure>`__ is a data structure consisting of a graph-theoretic tree with additional structure. CS trees are finite and have a `root node <http://en.wikipedia.org/wiki/Tree_(graph_theory)#rooted_tree>`_. As graphs, they are directed with edges pointing away from the root. Additionally, CS trees are ordered, meaning that every node has some specified order for its children.

This module works with trees of the latter kind. From this point onward, we use CS terminology (nodes, branches) rather than that of graph theory (vertex, edge or arc or arrow), and use the word *tree* to mean a CS tree.

.. testsetup:: 
	
	from thompson.trees import *
"""

from .drawing import *
from .constants import *

from collections import defaultdict, namedtuple
from fractions import Fraction

import svgwrite

__all__ = ['BinaryTree', 'DrawableTree']

class BinaryTree:
	"""The tree structure we require to work with elements of Thompsons's groups is implemented by :py:class:`BinaryTree`. A binary tree is a tree in which each node has at most 2 children.
	
	.. note::
		
		Whenever a node is expected, we use the value ``None`` to indicate that there is no such node. For example, if ``x`` is a node with no parent, then ``x.parent is None``.
	
	:ivar parent:	``None``: This node's parent.
	:ivar left:		``None``: This node's left child.
	:ivar right:	``None``: This node's right child."""
	
	def __init__(self, parent=None, right_child=False):
		"""Creates a new tree node. Attaches the new node to the `parent` if one is given. Otherwise, the new node is the root of a new tree.
		
		:param node parent:				The node to attach this node to, or `None` if this node is to have no parent.
		:param bool right_child: 		If True, attaches self as the right child of `parent`. If false, attaches self as the left child."""
		
		self.parent = parent
		self.left = None
		self.right = None
		
		if parent is not None:
			if right_child:
				parent.right = self
			else:
				parent.left = self
	
	def add_child(self, right_child=False):
		"""Creates a new node and adds it as a child of the current node.
		
		:param bool right_child:	If True, the new node is attached to the right; otherwise to the left.
		:return:					The newly created child node."""
		child = type(self)(parent=self, right_child=right_child)
		return child
	
	def is_root(self):
		"""Returns ``True`` if this node has no parent; otherwise ``False``."""
		return self.parent is None
	
	def is_leaf(self):
		"""Returns ``True`` if this node has no children; otherwise ``False``."""
		return self.left is None and self.right is None
	
	def num_leaves(self):
		"""Returns the number of leaves below this node.
		
		:rtype: int"""
		return sum(child.is_leaf() for child in self.walk())
	
	def num_children(self):
		"""Returns the number of children belonging to this node.
		
		:rtype: int"""
		return (self.left is not None) + (self.right is not None)
	
	def __iter__(self):
		"""Iterating over a node yields its left and right children. Missing children (``None``) are not included in the iteration."""
		if self.left is not None:
			yield self.left
		if self.right is not None:
			yield self.right
	
	def walk(self):
		"""The walk methods return generators that yield this node and its descendants, eventually iterating over the whole tree. The `three ways of traversing a tree depth-first <http://en.wikipedia.org/wiki/Tree_traversal#Types>`_ are available. If we don't care about the order of iteration, :py:meth:`walk` is a convenient shorthand for :py:meth:`walk_preorder`.
		
		.. todo::
			
			Currently implemented using recursion---may be more efficient using a stack/queue.
		"""
		#Normally I would just write walk = walk_preorder in the class body, but this lets me use Sphinx to document walk before walk_preorder.
		yield from self.walk_preorder()
		
	def walk_preorder(self):
		"""Yields the current node, then its left subtree, and then its right subtree."""
		yield self
		if self.left is not None:
			yield from self.left.walk_preorder()
		if self.right is not None:
			yield from self.right.walk_preorder()
	
	walk = walk_preorder
	
	def walk_in_order(self):
		"""Yields the current node's left subtree, then the current node itself, then the right subtree."""
		if self.left is not None:
			yield from self.left.walk_in_order()
		yield self
		if self.right is not None:
			yield from self.right.walk_in_order()
	
	def walk_postorder(self):
		"""Yields the left subtree, then the right subtree, and then finally the current node"""
		if self.left is not None:
			yield from self.left.walk_postorder()
		if self.right is not None:
			yield from self.right.walk_postorder()
		yield self
	
	@classmethod
	def from_string(cls, pattern):
		"""Creates a strict binary tree from a string (or any iterable) specifying the pattern of branches and leaves in the tree. A ``"1"`` represents a parent with two children, and a ``"0"`` represents a leaf. All other characters are ignored.
		
		The pattern specifies the tree in :py:meth:`pre-order <walk_preorder>` (since we create parents before either of their children).
		For example:
		
		- ``"100"``     denotes a single caret (``^``);
		- ``"1100100"`` denotes two carets with a common parent;
		- ``"11000"``   denotes a caret with another caret attached to its left.
		
		.. todo::
			
			Add pictures here."""
		#create the root
		pattern = cls.check_split_pattern(pattern)
		current = root = cls()
		for char in pattern:
			if char == "1":
				current = current.add_child()
			
			elif char == "0":
				while True: #poor man's do-while
					if current is root:
						assert root.num_children() != 1, "Root has exactly 1 child (expected 0 or 2 children). Pattern was %s." % pattern
						return root
					current = current.parent
					
					if current.left is not None and current.right is None:
						current = current.add_child(right_child=True)
						break
		return root
	
	"""TODO: ."""
	
	@staticmethod
	def check_split_pattern(pattern):
		"""Removes any characters from *pattern* that aren't ``'0'`` or ``'1'``. A few sanity checks are performed, to identify potential errors in the pattern. The method ensures that: 
		- the pattern describes an odd number of nodes,
		- the pattern does not start with a zero unless it describes a single node, and
		- the pattern has a sensible number of zeros and ones.
		
			>>> BinaryTree.check_split_pattern('12020')
			'100'
			>>> BinaryTree.check_split_pattern('1000')
			Traceback (most recent call last):
				...
			ValueError: Split patterns should be of odd length. Received '1000' of length 4.
			>>> BinaryTree.check_split_pattern('010')
			Traceback (most recent call last):
				...
			ValueError: Pattern '010' starts with a leaf but describes more than one node.
			>>> BinaryTree.check_split_pattern('10000')
			Traceback (most recent call last):
				...
			ValueError: Incorrect number of zeroes and ones. Received '10000' with 1 ones (expected 2).
		"""
		pattern = "".join(x for x in pattern if x in "01")
		
		if pattern[0] == '0' and len(pattern) > 1:
			raise ValueError("Pattern %r starts with a leaf but describes more than one node." % pattern)
		
		n = len(pattern)
		if n % 2 == 0:
			raise ValueError("Split patterns should be of odd length. Received %r of length %i." % (pattern, n))
		
		count = pattern.count("1")
		expected = (n - 1) // 2
		if count != expected:
			raise ValueError("Incorrect number of zeroes and ones. Received %r with %i ones (expected %i)." % (
				pattern, count, expected))
		return pattern
	
	def is_strict(self):
		"""Confusingly, different people use 'binary tree' to meant different things. We are interested in binary trees whose nodes:
		
		- have an order, so we can distinguish between left and right children, and
		- have exactly 0 or 2 children.
		
		Let us call such trees *strict* binary trees. This methods returns True iff the tree below this node is strictly binary.
		
			>>> BinaryTree.from_string("0").is_strict()
			True
			>>> BinaryTree.from_string("100").is_strict()
			True
			>>> x = BinaryTree()
			>>> left = x.add_child()
			>>> x.is_strict()
			False
		
		"""
		return all((child.left is None) == (child.right is None) for child in self.walk())
	
	def to_partition(self, partition=None, after=0, before=1):
		"""Strict binary trees represent a partition of an interval into dyadic sub-intervals. This method gives such a partition of :math:`[0, 1]`.
		
		.. todo::
			
			Use a stack/queue rather than recursion. Test with fraction classes."""
		if partition is None:
			partition = [Fraction(after), Fraction(before)]
		if self.is_leaf():
			return partition, after, before
		else:
			mid = (partition[after] + partition[after + 1]) / 2 
			partition.insert(after + 1, mid)
			partition, after, before = self.left.to_partition(partition, after, before+1)
			partition, after, before = self.right.to_partition(partition, after+1, before)
		
		return partition, after, before

Bounds = namedtuple('Bounds', 'min_x max_x height')

class DrawableTree(BinaryTree):
	"""The problem of drawing a binary tree (quite interesting in its own right) in the plane requires extra information than the basic tree structure provides (x and y coordinates, for instance). This subclass of :py:meth:`BinaryTree` provides methods which use that extra information to draw the tree tidily in the plane.
	
	The implementation of the :py:meth:`layout` algorithm was adapted from Bill Mill's article [Mill]_. Bill has made his code `available under the WTFPL <https://github.com/llimllib/pymag-trees>`_.
	
	:ivar x:		``0``			: This node's parent.
	:ivar y:		``0``			: This node's left child.
	:ivar bounds:	``(0, 0, 0)``	: A 3-tuple specifing the size of the subtree below and including this node. Entries are (*min_x*, *max_x*, *height*), storing the smallest/largest x-coordinate of any descendants and also the height of this node above its descendants.
	"""
	def __init__(self, parent=None, right_child=False):
		"""The same as :py:meth:`BinaryTree.__init__`. Some extra attributes needed for positioning are initialised."""
		super().__init__(parent, right_child)
		
		self.x = 0
		self.y = 0
		self.bounds = Bounds(0, 0, 0)
		
		self._offset = 0				#Used by layout methods to efficiently reposition a tree and its children.
		self._thread = None				#Used by layout_rt to determine when subtrees overlap.
	
	def layout(self):
		"""Uses Reingold and Tilfold's algorithm [RT]_ to calculate :math:`(x, y)` coordinates for each of the nodes. Among other aesthetic benefits, this algorithm always assigns integer coordinates to each node.
		"""
		self._setup_rt()
		self._add_offset()
		self._calc_bounds()
	
	def _setup_rt(self, depth=0):
		#1. First, the tree is traversed from the bottom up.
		for child in self:
			child._setup_rt(depth+1)
		
		self.y = depth
		#2. Leave leaves alone. Their parents will consider them a subtree and fix their position later.
		if self.is_leaf():
			self.x = 0
		else:
			assert self.num_children() == 2, "Expected exactly 0 or 2 children."
			self._fix_subtrees()
	
	def _fix_subtrees(self):
		#3. Otherwise, we're a parent. We have to move the right subtree enough to the right to ensure that it doesn't overlap the left subtree. To determine the distance, we calculate contours. A tree's left contour is the list of its leftmost nodes' x-coordinates at each depth; similarly for a right contour.
		
		#The contour function examines these contours and returns `sep' as before. Extra data is returned so we can set up threads for contour to function higher up the tree.
		
		li, ri, separation, loffset, roffset, lo, ro = _contour(self.left, self.right)
		#4. We have to move the trees separation+1 units apart to make sure no nodes share an x-coordinate. Once we've moved the tree to the right, the quantity in brackets will be the distance between self.left and self.right. This should be even, so that the parent (self) has an integer x coordinate.
		separation += 1
		if (self.right.x + separation - self.left.x) % 2 != 0:
			separation += 1
		
		self.right._offset = separation
		self.right.x += separation
		
		#5. We leave threads on any leaves on outer chords, so that chords can easily be retraced from higher in the tree. I'm not sure what the offsets are doing exactly other than allowing the algorithm to continue working futher up the tree.
		if not self.right.is_leaf():
			roffset += separation
		
		if ri and not li:							#For example the tree 10100 needs a thread on its left chord
			lo._thread = ri
			lo._offset = roffset - loffset
		elif li and not ri:							#For example the tree 11000 needs a thread on its right chord
			ro._thread = li
			ro._offset = loffset - roffset
		
		assert (self.left.x + self.right.x) % 2 == 0,  "Parent placed at non-integer coordinate"
		self.x =  (self.left.x + self.right.x) // 2
	
	def _next_left(self):
		return self._thread or self.left
	
	def _next_right(self):
		return self._thread or self.right 
	
	def _add_offset(self, offset_sum=0):
		#As we move up the tree in :py:meth:`layout`, we may find that we have to reposition subtrees. The repositioning information is stored in self._offset; this method uses that information."""
		self.x += offset_sum
		for child in self:
			child._add_offset(offset_sum + self._offset)
	
	def _calc_bounds(self):
		"""In one final post-order traversal of the tree, we compute the bounds attribute, which describes the position and height of the tree including and below the current node."""
		for child in self:
			child._calc_bounds()
		
		if self.is_leaf():
			min_x, max_x, height = self.x, self.x, 0
		else:
			#leaves have bounds = (0, 0, 0) already thanks to __init__
			min_x  = min(child.bounds.min_x  for child in self)
			max_x  = max(child.bounds.max_x  for child in self)
			height = max(child.bounds.height for child in self) + 1
		self.bounds = Bounds(min_x, max_x, height)
	
	def render(self, filename=None, leaf_names=None):
		"""Returns an SVG :py:class:`Group <svgwrite:svgwrite.container.Group>` whose contents are the drawing of this node and its descendants. 
		
		:param str filename:	If omitted, returns the group object. If specified, the group is added to a :class:`Drawing <svgwrite:svgwrite.drawing.Drawing>`. The drawing is then saved to *filename* and returned.
		:param leaf_names:		An optional list of names to be given to the leaves below this node.
		:type leaf_names: list of str
		"""
		if filename is not None:
			dwg, canvas = new_drawing(filename)
		g = svgwrite.container.Group()
		#The +1 is the extra space required to draw circles and not points.
		g.size = self.bounds.max_x - self.bounds.min_x + 1, self.bounds.height + 1
		
		#1. Draw the branches of the tree.
		for child in self.walk():
			if child.parent is None:
				continue
			start = Coord(child.parent)
			end = Coord(child)
			line = svgwrite.shapes.Line(start, end)
			g.add(line)
		
		#2. Draw all the nodes as circles
		i = 0
		for child in self.walk():
			#The order shouldn't matter---see the wikipedia page for tree traversal: leaves are always visited in the same order, because we visit the left subtree before the right.
			if leaf_names is not None and child.is_leaf():
				name = leaf_names[i]
				i += 1
			else:
				name = None
			g.add(child.render_node(name))
		
		#3. Offset the drawing so it aligns nicely with the grid.
		offset_group(g)
		
		if filename is not None:
			canvas.add(g)
			dwg.save()
			return dwg
		else:
			return g
	
	def render_node(self, name=None):
		"""Creates an SVG :py:class:`Circle <svgwrite:svgwrite.shapes.Circle>` representing this node.
		
		If *name* is given, returns a :py:class:`Group <svgwrite:svgwrite.container.Group>` containing the circle and
		*name* rendered as :py:class:`Text <svgwrite:svgwrite.text.Text>`. Otherwise, the circle is returned."""
		center = Coord(self)
		if name is None:
			circle = svgwrite.shapes.Circle(center, NODE_RADIUS)
		else:
			g = svgwrite.container.Group()
			circle = svgwrite.shapes.Circle((0, 0), NODE_RADIUS)
			g.translate(center)
			g.add(circle)
			
			text = svgwrite.text.Text(name, class_='centered')
			g.add(text)
		
		classes = []
		if self.is_leaf():
			classes.append("leaf")
		if self.parent and self is self.parent.left:
			classes.append("left")
		if classes:
			circle['class'] = " ".join(classes)
		
		if name is None:
			return circle
		else:
			return g

def _contour(left, right, max_sep=None, loffset=0, roffset=0, left_outer=None, right_outer=None):
	#See the comments for DrawableTree._fix_subtrees 
	#1. Compute the separation between the root nodes `left` and `right`, accounting for offsets.
	separation = left.x + loffset - (right.x + roffset)
	if max_sep is None or separation > max_sep:
		max_sep = separation
	
	if left_outer  is None: left_outer  = left
	if right_outer is None: right_outer = right
	
	lo = left_outer._next_left()					#Tracks the left contour of the left tree
	li = left._next_right()							#Tracks the right contour of the left tree
	ri = right._next_left()							#Tracks the left contour of the right tree
	ro = right_outer._next_right()					#Tracks the right contour of the right tree
	
	#2. If the inner chords continue, accumulate the offset from this depth and continue.
	if li and ri:
		loffset += left._offset
		roffset += right._offset
		return _contour(li, ri, max_sep, loffset, roffset, lo, ro)
		
	#3. If one of the trees has ended before the other, we can just return.
	return li, ri, max_sep, loffset, roffset, left_outer, right_outer