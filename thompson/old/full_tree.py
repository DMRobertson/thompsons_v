"""In graph theory, the word `tree <http://en.wikipedia.org/wiki/Tree_(graph_theory)>`__ means a connected graph with no cycles. In computer science (CS), a `tree <http://en.wikipedia.org/wiki/Tree_structure>`__ is a data structure consisting of a graph-theoretic tree with additional structure. CS trees are finite and have a `root node <http://en.wikipedia.org/wiki/Tree_(graph_theory)#rooted_tree>`_. As graphs, they are directed with edges pointing away from the root. Additionally, CS trees are ordered, meaning that every node has some specified order for its children.

This module works with trees of the latter kind. From this point onward we use the word *tree* to mean a CS tree.

.. testsetup:: 
	
	from thompson.full_tree import *
"""

class FullTree:
	"""A *k-ary* tree is a (CS) tree with the additional property that every node has at most *k* children. We call the parameter *k* the *arity* of the tree. Such a tree is called `full <http://en.wikipedia.org/wiki/K-ary_tree#Types_of_k-ary_trees>`_ if every node has either 0 or *k* leaves.
	
	This class implements full *k*-ary trees.
	
	:ivar arity:		The value *k*.
	:ivar parent:		This node's parent, or None if this node has no parent.
	:ivar children:		A list of this node's children. If the node is a leaf, its children is a list of *k* copies of ``None``.
	:ivar data:			Nodes can carry a reference to any Python object they might stand for. Initially None.
	"""
	# __slots__ = ("arity", "parent", "children", "data")
	
	#Creation
	def __init__(self, arity):
		"""Creates a new node in a full *k*-ary tree. After creation, the node is a leaf."""
		self.arity = arity
		self.parent = None
		self.children = [None]*arity
		self.data = None
	
	#Tests
	def is_root(self):
		"""Returns True if this node has no parent, otherwise False.
		
			>>> node = FullTree(2)
			>>> node.is_root()
			True
			>>> node.expand()
			>>> node.children[0].is_root()
			False
		"""
		return self.parent is None
	
	def is_leaf(self):
		"""Returns True if this node has no children, otherwise False.
		
			>>> node = FullTree(2)
			>>> node.is_leaf()
			True
			>>> node.expand()
			>>> node.is_leaf(), node.children[0].is_leaf()
			(False, True)
		"""
		#Assumption: the tree is full.
		return self.children[0] is None
	
	def __len__(self):
		"""The length of a node is 0 if it is a leaf, or else the arity *k*. 
		
			>>> node = FullTree(2)
			>>> len(node)
			0
			>>> node.expand()
			>>> len(node)
			2
		
		"""
		#Assumption: the tree is full.
		return 0 if self.is_leaf() else self.arity
	
	def __bool__(self):
		#If we do not define this method, Python will assume bool(x) == len(x) > 0. But it is useful to have any node be truthy compared to None which is falsey. (e.g. for short circuiting).
		return True
	
	#Traversal
	def __iter__(self):
		"""Iterating over a node yields nothing if it is a leaf; otherwise iterating yields the node's children."""
		if not self.is_leaf():
			yield from self.children
	
	def walk(self):
		"""The walk methods return generators that yield this node and its descendants, eventually iterating over the whole tree. There are `three ways of traversing a binary (2-ary) tree depth-first <http://en.wikipedia.org/wiki/Tree_traversal#Types>`_ . Both pre- and post-order traversal generalise to k-ary trees.
		
		If we don't care about the order of iteration, :py:meth:`walk` is a convenient shorthand for :py:meth:`walk_preorder`.
		
			>>> root = FullTree(3); root.expand()
			>>> for child in root: child.expand()
			>>> sum(1 for child in root.walk())
			13
		"""
		#Normally I would just write walk = walk_preorder in the class body, but this lets me use Sphinx to document walk before walk_preorder.
		yield from self.walk_preorder()
	
	def walk_preorder(self):
		"""Yields the current node and then its children.
		
			>>> root = FullTree(3); root.expand()
			>>> node = next(root.walk_preorder())
			>>> node.is_root(), node.is_leaf()
			(True, False)
		"""
		yield self
		for child in self:
			yield from child.walk_preorder()
	
	def walk_postorder(self):
		"""Yields the current node's children before itself.
		
			>>> root = FullTree(3); root.expand()
			>>> node = next(root.walk_postorder())
			>>> node.is_root(), node.is_leaf()
			(False, True)
		"""
		for child in self:
			yield from child.walk_postorder()
		yield self
	
	def __str__(self):
		"""Returns a string representing the pattern of branches and leaves in this tree in pre-order traversal. ``'1'`` represents a branch and ``'0'`` represents a leaf.
		
			>>> x = FullTree(3); str(x)
			'0'
			>>> x.expand(); str(x)
			'1000'
			>>> x.children[2].expand(); str(x)
			'1001000'
			>>> x.children[0].expand(); str(x)
			'1100001000'
		"""
		if self.is_leaf():
			return '0'
		return '1' + "".join(child.__str__() for child in self)
	
	#Modifications
	def _add_child(self, index):
		r"""Creates a new tree node and adds it to position *index* in this node's list of children.
		
		:raises IndexError: if *index* is outside of the range :math:`0, \dotsc, \text{arity} - 1`.
		"""
		child = type(self)(self.arity)
		child.parent = self
		self.children[index] = child
		return child
	
	def expand(self):
		"""If called on a leaf, turns this leaf into a branch with *k* newly-created children.
		
		:raises ValueError: if called on a node which is not a leaf.
		"""
		if not self.is_leaf():
			raise ValueError("This node has already been expanded.")
		for i in range(self.arity):
			self._add_child(i)
	
	def replace_with_child(self, index):
		r"""Let *S* (for self) stand for the current node. This method:
		
		- removes the child node *C* with index *index* from *S* ;
		- removes all other children from *S*.
		- if *S* has a parent *P* :
		
			- replaces *S* with *C* in *P*'s list of children;
			- sets the parent of *C* to be *P*.
		- Finally, the method returns *C*.
		
		:raises IndexError: if *index* is outside of the range :math:`0, \dotsc, \text{arity} - 1`.
		:raises ValueError: if *S* is a leaf.
		:returns: the child node *C*.
		"""
		replacement = self.children[index]
		if replacement is None:
			raise ValueError("Cannot call replace_with_child on a leaf.")
		for child in self:
			if child is replacement:
				continue
			child.discard()
		self.children = [None] * self.arity
		
		parent = replacement.parent = self.parent
		self.parent = None
		if parent is not None:
			j = parent.children.index(self)
			parent.children[j] = replacement
		return replacement
	
	def discard(self):
		"""Removes the reference
		- from this node to its parent,
		- from any descendant of this node to their parents,
		- from any descendant of this node to their children, and
		- from this node to its children.
		
		This isn't strictly necessary but it's good bookkeeping. See the note about :meth:`__del__ <py3:object.__del__>` in the Python documentation.
		
		*NB:* this method does NOT remove the reference from ``self.parent`` to ``self`` in ``self.parent.children``.
		"""
		self.parent = None
		for child in self:
			child.discard()
		self.children = [None] * self.arity
