class Word(list):
	def letter(self, index=1):
		r"""Writes the letter :math:`x_\text{index}` at the end of this word.
		
		:raises ValueError: if index is not in the range 1, ..., *alphabet_size*."""
		if not 1 <= index <= self.alphabet_size:
			raise ValueError("Index {} is not in the range 1 to alphabet size {}".format(index, self.alphabet_size))
		self.append(index)
		return self
	
	# x = letter
	
	def descend(self, index, power=1):
		r"""Writes the descend operator :math:`\alpha_\text{index}` at the end of this word *power* times.
		
		:raises ValueError: if index is not in the range 1, ..., *arity*."""
		if not 1 <= index <= self.arity:
			raise ValueError("Index {} is in the range 1 to arity {}".format(index, self.arity))
		for _ in range(power):
			self.append(-index)
		return self
	
	# alpha = descend
	
	def contract(self):
		r"""Writes the contraction operator :math:`\lambda` at the end of this word."""
		self.append(0)
		return self
	
	# lambda_ = contract
	
	@classmethod
	def initial_segment(cls, u, v):
		r"""Returns True if :math:`u` is an initial segment of :math:`v` or vice versa; otherwise returns False.
		
		Def 3.15: :math:`u` is an initial segment of :math:`v` if :math:`v = u\Gamma`, where :math:`\Gamma` is some string of alphas. Call :math:`u` a *proper* initial segment of :math:`v` if :math:`\Gamma` is not the empty string, *i.e.* if :math:`u \neq v`.
		
			>>> #completely different words
			>>> u = Word.from_string("x1 a2 a1 a1 a2")
			>>> v = Word.from_string("x2 a1 a2 a1 a2")
			>>> Word.initial_segment(u, v)
			False
			>>> #v starts with u
			>>> u = Word.from_string("x a1 a1 x2 a1 a2 L a1 a2")
			>>> v = Word.from_string("x a1 a1 a2 a2 a1")
			>>> Word.initial_segment(u, v)
			True
			>>> #v starts with u but has a lambda contraction later
			>>> u = Word.from_string("x a1 a2")
			>>> v = Word.from_string("x a1 a2 x a2 a2 L")
			>>> Word.initial_segment(u, v)
			False
		"""
		if len(u) > len(v):
			u, v = v, u
		u.standardise()
		v.standardise()
		return all(a == b for a, b in zip(u, v)) and all(b < 0 for b in v[len(u):])
		
			