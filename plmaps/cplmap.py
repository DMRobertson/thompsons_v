from  math     import ceil, floor
from fractions import Fraction

from .plmap import PLMap, PL2
from .util  import grad, lerp

class CPLMap(PLMap):
	r"""
	A :class:`PLMap` specialised to represent a PL map of the circle :math:`[0, 1] / \{ 0 \sim 1 \}`.
	Such an element must be input as if it were a PLMap :math:`[a, a+1] \to [b, b+1]`.
	
	"""
	
	def __init__(self, domain, range):
		domain = [Fraction(d) for d in domain]
		range  = [Fraction(r) for r in range ]
		super()._validate_inputs(domain, range)
		for name, list in (("domain", domain), ("range", range)):
			if list[-1] - list[0] != 1:
				raise ValueError(name + " doesn't have width 1")
		
		self.domain, self.range, self.gradients = self._normalise_breakpoints(domain, range)
	
	def _normalise_breakpoints(self, domain, range):
		#Between which indices does x = 0 lie?
		for i, _ in enumerate(domain):
			if domain[i - 1] > domain[i]:
				break
		exact = domain[i] == 0
		
		#normalise x coordinates
		for j, _ in enumerate(domain):
			domain[j] %= 1
		if exact:
			domain[i - 1] = Fraction(1)
		
		#normalise y coordinates
		offset = floor(range[i])
		for j, _ in enumerate(range):
			range[j] -= offset
		
		#add in a (0, y) and later a (1, y+1) coordinate if neccesary
		if exact:
			new_domain = []
			new_range  = []
		else:
			#must be in between x=domain[i -1] and x=domain[i]
			#where does the line segment equal one?
			initial_y = lerp(
				ceil(domain[i - 1]),
				domain[i - 1], domain[i] + 1,
				range [i - 1], range [i]
			)
			new_domain = [0]
			new_range  = [initial_y]
		
		new_domain.extend(domain[i:])
		new_range.extend(range[i:])
		new_domain.extend(domain[:i])
		new_range.extend(r + 1 for r in range[:i])
		
		if not exact:	
			new_domain.append(Fraction(1))
			new_range.append(initial_y + 1)
		return super()._normalise_breakpoints(new_domain, new_range)
	
	def _validate_multiplication(self, other):
		pass
	
	def _compute_product_breakpoints(self, other):
		#Slightly awkward hack here...
		domain, range = super()._compute_product_breakpoints(other)
		if range[-1] == range[0]:
			range[-1] += 1
		return domain, range
	
	def image(self, x):
		"""Take in anything; spew out something in [0, 1)"""
		return super().image(x % 1) % 1
	
	def inverse_image(self, y):
		"""Take in anything; spew out something in [0, 1)"""
		#shift by an integer so that y + n \in [r, r+1]
		r = self.range[0]
		offset = (y - r) % 1
		return super().inverse_image(r + offset)
		
class CPL2(CPLMap, PL2):
	pass