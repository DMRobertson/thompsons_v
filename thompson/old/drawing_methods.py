	def layout_knuth(self, depth=0, next=None):
		"""Assigns each node to its own column."""
		self.y = depth
		
		"""next[0] tracks the next available x coordinate. We store this quantity in a list so that changes to it are recognised further up the tree."""
		if next is None:
			next = [0]
		
		if self.left:
			self.left.layout_knuth(depth+1, next)
		
		self.x = next[0]
		next[0] += 1
		
		if self.right:
			self.right.layout_knuth(depth+1, next)
	
	def layout_ws1(self, depth=0, nexts=None):
		"""Wetherell and Shannon's first strategy. Like Knuth's layout, but 'next' is tracked for each row, so that we can use the minimal amount of horizontal space."""
		if nexts is None:
			nexts = defaultdict(int)
		
		self.x = nexts[depth]
		nexts[depth] += 1
		self.y = depth
		
		for child in self:
			child.layout_ws1(depth+1, nexts)
	
	def layout_ws2(self):
		"""Wetherell and Shannon's second strategy. Centres parents above their children and ensures correct spacing."""
		self.setup_ws2()
		self.add_offset()
	
	def setup_ws2(self, depth=0, nexts=None, offsets=None):
		self.y = depth
		if nexts   is None: nexts   = defaultdict(int)			#next available position
		if offsets is None: offsets = defaultdict(int)			#sum of offsets applied at this depth
		
		"""1. First, the tree is traversed from the bottom up."""
		for child in self:										
			child.setup_ws2(depth+1, nexts, offsets)
		
		if self.is_leaf():
			self.x = nexts[depth]								#Just put leaves in the next available spot.
		else:
			assert self.num_children() == 2, "Expected a non-leaf node to have exactly 2 children."
			"""2. This means that we can centre a parent above its children by taking the average of the children's x-coordinates."""
			self.x = (self.left.x + self.right.x)/2
		
		"""3. Sometimes centring a node moves it left of the next available x-coordinate. (E.g. the tree 10100.) To compensate, we note that it and its descendants should be moved the same distance to the right in a second pass. """
		if offsets[depth] + self.x < nexts[depth]:				#if the current position is too far left
			offsets[depth] = nexts[depth] - self.x				#increase offset to position parent at nexts[depth]
		
		"""TODO: What does this if statement do?"""
		if not self.is_leaf():
			self.x += offsets[depth]
		
		"""4. Equally parents are sometimes centred to the right of the next available position (E.g. the tree 11000). To compensate, we must increase the appropriate value in 'nexts' to indicate that space has been taken up."""
		nexts[depth] = self.x + 2
		self.offset = offsets[depth]
	
	def fix_subtrees_naive(self):
		"""3. Otherwise, we're a parent. We have to move the right subtree enough to the right to ensure that it doesn't overlap the left subtree. To determine the distance, we calculate contours. A tree's left contour is the list of its leftmost nodes' x-coordinates at each depth; similarly for a right contour."""
		lr = self.left.contour_naive(right_edge=True)		#left tree's right contour
		rl = self.right.contour_naive(right_edge=False)		#right tree's left contour
		
		"""At each depth, compute the distance between the left subtree's right contour and the right subtree's left contour. Call the largest of these values `separation'. We have to move the trees separation+1 units apart to make sure no nodes share an x-coordinate."""
		separation = max(x-y for x,y in zip(lr, rl)) + 1
		"""Once we've moved the tree to the right, the quantity in brackets will be the distance between self.left and self.right. This should be even, so that the parent (self) has an integer x coordinate."""
		if (self.right.x + separation - self.left.x) % 2 != 0:
			separation += 1
		
		"""4. Add the separation and leave a note to offset any descendants of the right tree."""
		self.right.x += separation
		self.right.offset += separation
		
		assert (self.left.x + self.right.x) % 2 == 0,  "Parent placed at non-integer coordinate"
		self.x = (self.left.x + self.right.x) // 2
	
	def contour_naive(self, right_edge=False, depth=0, cont=None):
		"""A top-down method to compute contours of a subtree. We store the smallest (resp. largest) x-coordinates seen at each depth for a left (resp. right) contour. This potentially checks every node in the map; TR's more efficient algorithm involves `threads'."""
		if cont is None:
			cont = [self.x]
		elif len(cont) < depth + 1:
			cont.append(self.x)
		elif right_edge and cont[depth] < self.x:
			cont[depth] = self.x
		elif not right_edge and self.x < cont[depth]:
			cont[depth] = self.x
		
		for child in self:
			child.contour_naive(right_edge, depth+1, cont)
		
		return cont