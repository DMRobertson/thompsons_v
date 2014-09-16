"""In graph theory, the word `tree <http://en.wikipedia.org/wiki/Tree_(graph_theory)>`__ means a connected graph with no cycles. In computer science (CS), a `tree <http://en.wikipedia.org/wiki/Tree_structure>`__ is a data structure consisting of a graph-theoretic tree with additional structure. CS trees are finite and have a `root node <http://en.wikipedia.org/wiki/Tree_(graph_theory)#rooted_tree>`_. As graphs, they are directed with edges pointing away from the root. Additionally, CS trees are ordered, meaning that every node has some specified order for its children.

This module works with trees of the latter kind. From this point onward, we use CS terminology (nodes, branches) rather than that of graph theory (vertex, edge or arc or arrow), and use the word *tree* to mean a CS tree.

.. testsetup:: 
	
	from thompson.trees import *
"""

from .drawing import *
from .constants import *

from collections import defaultdict, namedtuple
from fractions import Fraction
import random

import svgwrite

__all__ = ['BinaryTree', 'DrawableTree', 'random_tree']

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
		
		if isinstance(parent, str):
			raise TypeError('Parent should be a tree or None, not a string. did you mean to use BinaryTree.from_string()?')
		
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
	
	#TODO merge these two methods
	
	def set_child(self, child, right_child=False):
		"""Sets the current node's child to be the given node. You should use this method rather than alter ``self.left`` or ``self.right`` directly, because ``child.parent`` needs to be updated too.
		
		Beware: cycles can be created with this method, meaning that the graph isn't a tree any more. Any method that traverses the tree will then attempt to loop infinitely! For instance:
		
			>>> node = BinaryTree()
			>>> node.set_child(node)
			>>> node.set_child(node, right_child=True)
			>>> for descendant in node.walk(): pass
			Traceback (most recent call last):
				...
			RuntimeError: maximum recursion depth exceeded
		"""
		if right_child:
			self.right = child
		else:
			self.left = child
		child.parent = self
	
	def detach_children(self, recursive=True):
		"""Detaches this node's children (if any) from itself. If *recursive* is True, the descendants of this node have their children detached too, `allowing the Python garbage collector <py3:object.__del__>` to remove them from memory. Otherwise the subtrees below this node are left intact."""
		if recursive:
			for child in self:
				child.detach_children(recursive=True)
		for child in self:
			child.parent = None
		self.left = self.right = None
	
	def is_root(self):
		"""Returns ``True`` if this node has no parent; otherwise ``False``.
		
			>>> x = BinaryTree.from_string("100")
			>>> x.is_root()
			True
			>>> x.left.is_root(), x.right.is_root()
			(False, False)"""
		return self.parent is None
	
	def is_leaf(self):
		"""Returns ``True`` if this node has no children; otherwise ``False``.
		
			>>> x = BinaryTree.from_string("100")
			>>> x.is_leaf()
			False
			>>> x.left.is_leaf(), x.right.is_leaf()
			(True, True)
		"""
		return self.left is None and self.right is None
	
	def is_trivial(self):
		"""Returns ``True`` if this node has no parent and no children; otherwise ``False``.
		
			>>> BinaryTree().is_trivial()
			True
			>>> BinaryTree.from_string("100").is_trivial()
			False
		"""
		return self.parent is None and self.left is None and self.right is None
	
	def is_left_child(self):
		"""Returns ``True`` if this node has a parent and if it is the left child of the parent.
		
			>>> x = BinaryTree()
			>>> x.is_left_child()
			False
			>>> x = BinaryTree.from_string("100")
			>>> x.left.is_left_child(), x.right.is_left_child()
			(True, False)
		"""
		return self.parent is not None and self is self.parent.left
	
	def is_right_child(self):
		"""Returns ``True`` if this node has a parent and if it is the right child of the parent.
		
			>>> x = BinaryTree()
			>>> x.is_right_child()
			False
			>>> x = BinaryTree.from_string("100")
			>>> x.left.is_right_child(), x.right.is_right_child()
			(False, True)
		"""
		return self.parent is not None and self is self.parent.right
	
	def __len__(self):
		"""We define the length of a tree node to be the number of non-``None`` children it possesses.
		
			>>> len(BinaryTree())
			0
			>>> len(BinaryTree.from_string("100"))
			2
			>>> root = BinaryTree()
			>>> child = root.add_child()
			>>> len(root)
			1
		"""
		return (self.left is not None) + (self.right is not None)
	
	def __bool__(self):
		#To evaluate bool(x), Python looks for x.__bool__(). Failing that, Python tries x.__len__() and returns len(x) > 0. This means that python would treat leaves as being False!
		#In some places like _next_left(), we use the short circuiting of or to avoid returning None where possible.
		#For this to work, non-None instances should be True.
		return True
	
	def num_leaves(self):
		"""Returns the number of leaves below and including this node.
		
			>>> BinaryTree().num_leaves() #just a single leaf
			1
			>>> BinaryTree.from_string("100").num_leaves()
			2
		"""
		return sum(child.is_leaf() for child in self.walk())
	
	def __eq__(self, other):
		"""Two trees are considered equal (isomorphic) if they have the same local structure at every node.
		
			>>> BinaryTree() == BinaryTree() #different instances, same structure
			True
			>>> BinaryTree.from_string("100") == BinaryTree.from_string("10100")
			False
			>>> BinaryTree.from_string("10100") != BinaryTree.from_string("11000") #mirror image
			True
			>>> BinaryTree() == "hello"
			False
		"""
		if not isinstance(other, BinaryTree):
			return NotImplemented
		
		#Ensure self.left exists iff other.left exists
		if (self.left is None) != (other.left is None):
			return False
		#If they exist, check the left subtrees are equal
		if (self.left is not None) and self.left != other.left:
			return False
		
		#Ensure self.right exists iff other.right exists
		if (self.right is None) != (other.right is None):
			return False
		#If they exists, check the left subtrees are equal
		if (self.right is not None) and self.right != other.right:
			return False
		return True
	
	def __ne__(self, other):
		if not isinstance(other, BinaryTree):
			return NotImplemented
		return not (self == other)
	
	def __iter__(self):
		"""Iterating over a node yields its left and right children. Missing children (``None``) are not included in the iteration."""
		if self.left is not None:
			yield self.left
		if self.right is not None:
			yield self.right
	
	def walk(self):
		"""The walk methods return generators that yield this node and its descendants, eventually iterating over the whole tree. The `three ways of traversing a tree depth-first <http://en.wikipedia.org/wiki/Tree_traversal#Types>`_ are available. If we don't care about the order of iteration, :py:meth:`walk` is a convenient shorthand for :py:meth:`walk_preorder`.
		"""
		#Currently implemented using recursion---may be more efficient using a stack/queue.
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
		
		The pattern specifies the tree in :meth:`pre-order <walk_preorder>` (since we create parents before either of their children). Some sanity checks are made by :meth:`check_split_pattern`.
		For example:
		
		- ``"100"``     denotes a single caret (``^``);
		- ``"1100100"`` denotes two carets with a common parent;
		- ``"11000"``   denotes a caret with another caret attached to its left.
		
		.. figure:: /examples/trees/BinaryTree_from_string.svg
		
			The trees formed by the strings above. [:download:`Source code </examples/trees/BinaryTree_from_string.py>`].
		
		:raises ValueError: if :meth:`check_split_pattern` finds an error in the pattern
		:raises ValueError: if the pattern is ill-formed and ends prematurely. For example, "10001" appears to describe a tree with two branches and three leaves. After reading "100", the algorithm will stop and the final two characters will be ignored. 
		
			>>> BinaryTree.from_string("10001")
			Traceback (most recent call last):
				...
			ValueError: Pattern 10001 ended prematurely after character #3.
		"""
		pattern = cls.check_split_pattern(pattern)
		break_outer = False
		#create the root
		current = root = cls()
		for i, char in enumerate(pattern):
			if char == "1":
				current = current.add_child()
			elif char == "0":
				while True: #poor man's do-while
					if current is root:
						assert len(root) != 1, "Root has exactly 1 child (expected 0 or 2 children). Pattern was %s." % pattern
						break_outer = True
						break
					current = current.parent
					
					if current.left is not None and current.right is None:
						current = current.add_child(right_child=True)
						break
			if break_outer:
				break
		
		if i + 1 != len(pattern):
			raise ValueError('Pattern {} ended prematurely after character #{}'.format(pattern, i+1))
		return root
	
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
		
			>>> partition, _, _ = BinaryTree.from_string("11000").to_partition()
			>>> partition
			[Fraction(0, 1), Fraction(1, 4), Fraction(1, 2), Fraction(1, 1)]
			>>> print(*(str(x) for x in partition), sep=', ')
			0, 1/4, 1/2, 1
			
		:return: a list of :class:`Fractions <py3:fractions.Fraction>` forming the partition.
		"""
		#todo. Use a stack/queue rather than recursion.
		#TODO. I have no idea why this works, so this could probably do with a rewrite.
		if partition is None:
		
			partition = [Fraction(after), Fraction(before)]
		if self.is_leaf():
			return partition, after, before
		else:
			mid = (partition[after] + partition[after + 1]) / 2 
			partition.insert(after + 1, mid)
			partition, after, before = self.left.to_partition(partition, after, before+1)
			partition, after, before = self.right.to_partition(partition, after+1, before)
		
		if after == 0 and before == 1:
			assert len(partition) == self.num_leaves(), "Partition {} size does not match the number of leaves ({}).".format(partition, self.num_leaves())
		return partition, after, before
	
	def leaves(self, perm=None):
		"""Returns a depth-first list of leaves below this node.
		
		If a permutation *perm* is given, it is applied before returning the list. The resulting list will be in label order rather than traversal order."""
		size = self.num_leaves()
		if perm is not None and size != perm.size:
			raise ValueError("{} is of size {}, but the tree has {} leaves.".format(repr(perm), perm.size, self.num_leaves()))
		
		leaves = [node for node in self.walk() if node.is_leaf()]
		if perm is not None:
			new_list = [None] * size
			for i, image in perm:
				new_list[i-1] = leaves[image-1]
			leaves = new_list
		return leaves

Bounds = namedtuple('Bounds', 'min_x max_x height')

class DrawableTree(BinaryTree):
	"""The problem of drawing a binary tree (quite interesting in its own right) in the plane requires extra information than the basic tree structure provides (x and y coordinates, for instance). This subclass of :py:meth:`BinaryTree` provides methods which use that extra information to draw the tree tidily in the plane.
	
	The implementation of the :py:meth:`layout` algorithm was adapted from Bill Mill's article [Mill]_. Bill has made his code `available under the WTFPL <https://github.com/llimllib/pymag-trees>`_.
	
	:ivar x:		``0``
	:ivar y:		``0``
	:ivar bounds:	``(0, 0, 0)``. A 3-tuple specifing the position and height of the subtree below and including the current node. Automatically determined when :meth:`layout()` is called.
	"""
	def __init__(self, parent=None, right_child=False):
		"""The same as :py:meth:`BinaryTree.__init__`. Some extra attributes needed for positioning are initialised."""
		super().__init__(parent, right_child)
		self._reset_drawing_info()
	
	def _reset_drawing_info(self):
		"""Resets all state used to draw the tree."""
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
		assert self.bounds.min_x == 0, "Leftmost tree coordinate is x = %i (expected x=0)" % self.bounds.min_x
	
	def _setup_rt(self, depth=0):
		"""Traverses the tree in post-order to assign x-coordinates to the tree's nodes. Note that :meth:`_add_offset` should be called after this method in order to apply the offsets."""
		#If case we have modified the tree since a previous layout() call.
		self._reset_drawing_info()
		
		#We cannot determine a node's x-coordinate until we know those of its children.
		for child in self:
			child._setup_rt(depth+1)
		
		#The y-coordinates are easy
		self.y = depth
		if self.is_leaf():
			self.x = 0
		elif len(self) == 1:
			child = self.left or self.right
			self.x = child.x
		else:
			self.x = self._fix_subtrees()
	
	def _fix_subtrees(self):
		#3. Otherwise, we're a parent. We have to move the right subtree enough to the right to ensure that it doesn't overlap the left subtree. To determine the distance, we calculate contours. A tree's left contour is the list of its leftmost nodes' x-coordinates at each depth; similarly for a right contour.
		
		#The contour function examines these contours and returns `sep' as before. Extra data is returned so we can set up threads for contour to function higher up the tree.
		
		li, ri, separation, loffset, roffset, lo, ro = _contour(self.left, self.right)
		#4. We have to move the trees separation+1 units apart to make sure no nodes share an x-coordinate. Once we've moved the tree to the right, the quantity in brackets will be the distance between self.left and self.right. This should be even, so that the parent (self) has an integer x coordinate.
		separation += 1
		if (self.right.x + separation - self.left.x) % 2 != 0:
			separation += 1
		if separation > 0:
			self.right._offset = separation
			self.right.x += separation
		#Handling this case separately ensures that the leftmost x-coordinate is zero.
		elif separation < 0:
			self.left._offset = -separation
			self.left.x -= separation
		#5. We leave threads on any leaves on outer chords, so that chords can easily be retraced from higher in the tree. I'm not sure what the offsets are doing exactly other than allowing the algorithm to continue working further up the tree.
		if not self.right.is_leaf():
			roffset += separation
		
		if ri is not None and li is None:				#For example the tree 10100 needs a thread on its left chord
			lo._thread = ri
			lo._offset = roffset - loffset
		elif li is not None and ri is None:				#For example the tree 11000 needs a thread on its right chord
			ro._thread = li
			ro._offset = loffset - roffset
		
		assert (self.left.x + self.right.x) % 2 == 0,  "Parent placed at non-integer coordinate"
		return (self.left.x + self.right.x) // 2
	
	def _next_left(self):
		return self._thread or self.left or self.right
	
	def _next_right(self):
		return self._thread or self.right or self.left
	
	def _add_offset(self, offset_sum=0):
		"""As we move up the tree in :py:meth:`layout`, we may find that we have to reposition subtrees so they don't overlap. The repositioning information is stored in self._offset; this method uses that information."""
		self.x += offset_sum
		for child in self:
			child._add_offset(offset_sum + self._offset)
	
	def _calc_bounds(self):
		"""In one final post-order traversal of the tree, we compute the bounds attribute. This is a 3-tuple specifing the position and height of the subtree below and including the current node.
		
		Entries are (*min_x*, *max_x*, *height*), storing the smallest/largest x-coordinate of any descendants and also the height of this node above its descendants. Entries can be accessed by index or attribute.
		"""
		for child in self:
			child._calc_bounds()
		
		if self.is_leaf():
			min_x, max_x, height = self.x, self.x, 0
		else:
			min_x  = min(child.bounds.min_x  for child in self)
			max_x  = max(child.bounds.max_x  for child in self)
			height = max(child.bounds.height for child in self) + 1
		self.bounds = Bounds(min_x, max_x, height)
	
	@creates_svg
	def render(self, leaf_names=None, **kwargs):
		"""Returns an SVG :py:class:`Group <svgwrite:svgwrite.container.Group>` whose contents are the drawing of this node and its descendants. 
		
		:param leaf_names: An optional list of names to be given to the leaves below this node. The leaves are specified depth-first, left to right---the same order as all the :meth:`BinaryTree.walk` methods.
		:type leaf_names:  list of str
		
		.. figure:: /examples/trees/DrawableTree_render.svg
		
			**Example**. A randomly generated tree rendered to an SVG. [:download:`Source code <examples/trees/DrawableTree_render.py>`]
		"""
		g = svgwrite.container.Group()
		
		#1. Draw the branches of the tree.
		for child in self.walk():
			if child is self:
				continue
			start = Coord(child.parent)
			end = Coord(child)
			line = svgwrite.shapes.Line(start, end)
			g.add(line)
			
			if child._thread:
				start = Coord(child)
				end = Coord(child._thread)
				line = svgwrite.shapes.Line(start, end, class_='thread debug')
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
		
		#3. Before finishing, offset the group so it aligns nicely with the grid.
		#Add 1 to width and height to give extra room for the node's radii. See the documentation of set_size for an example.
		size = Coord(self.bounds.max_x - self.bounds.min_x + 1, self.bounds.height + 1)
		set_size(g, size, offset=Coord(0.5, 0.5))
		return g
	
	def render_node(self, name=None):
		"""Creates an SVG :py:class:`Circle <svgwrite:svgwrite.shapes.Circle>` representing this node.
		
		If *name* is given, returns a :py:class:`Group <svgwrite:svgwrite.container.Group>` containing the circle and
		*name* rendered as :py:class:`Text <svgwrite:svgwrite.text.Text>`. Otherwise, the circle is returned."""
		
		#TODO Maybe render the nodes deeper into the tree with a smaller radius? Or scale the y-axis a little
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

def _contour(left, right, max_sep=float('-inf'), loffset=0, roffset=0, left_outer=None, right_outer=None):
	#See the comments for DrawableTree._fix_subtrees 
	#1. Compute the separation between the root nodes `left` and `right`, accounting for offsets.
	separation = left.x + loffset - (right.x + roffset)
	max_sep = max(max_sep, separation)
	
	if left_outer  is None: left_outer  = left
	if right_outer is None: right_outer = right
	
	lo = left_outer._next_left()					#Tracks the left contour of the left tree
	li = left._next_right()							#Tracks the right contour of the left tree
	ri = right._next_left()							#Tracks the left contour of the right tree
	ro = right_outer._next_right()					#Tracks the right contour of the right tree
	
	#2. If the inner chords continue, accumulate the offset from this depth and continue.
	if li is not None and ri is not None:
		loffset += left._offset
		roffset += right._offset
		return _contour(li, ri, max_sep, loffset, roffset, lo, ro)
		
	#3. If one of the trees has ended before the other, we can just return.
	return li, ri, max_sep, loffset, roffset, left_outer, right_outer
def random_tree(num_leaves=None):
	"""Returns a :class:`DrawableTree` with *num_leaves* leaves. If *num_leaves* is omitted, it defaults to a random integer in the range 5, ..., 15.
	
	The pattern of branches is randomly determined.
	
	I first saw this implemented in Roman Kogan's `nvTrees applet <http://www.math.tamu.edu/~romwell/nvTrees/>`_.
	
	.. doctest::
		
		>>> from random import randint
		>>> x = randint(1, 20)
		>>> random_tree(x).num_leaves() == x
		True
	"""
	if num_leaves is None:
		num_leaves = random.randint(5, 15)
	
	subtrees = ["0"] * num_leaves
	indices = list(range(num_leaves))
	
	for _ in range(num_leaves-1):
		#Pick two subtrees at random
		i, j = random.sample(indices, 2)
		#and join them together
		subtrees[i] = "1" + subtrees[i] + subtrees[j]
		del indices[-1], subtrees[j]
	
	return DrawableTree.from_string(subtrees[0])