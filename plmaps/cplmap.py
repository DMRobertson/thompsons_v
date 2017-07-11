import builtins
from math      import ceil, floor
from fractions import Fraction

from .plmap import PLMap, PL2, linear_superclass
from .util  import ends, grad, lerp

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
	
	@staticmethod
	def _translate_near_zero(list):
		#Assume list is an increasing sequence of ints with list[0] + 1 == list[-1].
		#Modify the list in-place, so that each entry lies in [0, 2).
		#With the same remainder mod 1
		m = floor(min(list))
		for i, _ in enumerate(list):
			list[i] -= m
	
	def _normalise_breakpoints(self, domain, range):
		#_translate_near_zero modifies in place
		domain = list(domain)
		self._translate_near_zero(domain)
		#Between which indices does x = 1 lie?
		for i in builtins.range(1, len(domain)):
			if domain[i-1] < 1 <= domain[i] :
				break
		exact = domain[i] == 1
		
		new_domain = []
		new_range  = []
		
		image_of_one = lerp(
			1,
			domain[i - 1], domain[i],
			range [i - 1], range [i]
		)
		if not exact:
			#Add in (0, image of 0)
			new_domain.append(Fraction(0))
			new_range .append(image_of_one - 1)
		#Add in all values greater than one, after reducing them mod 1
		new_domain.extend( d % 1 for d in domain[i:] )
		new_range .extend( r - 1 for r in range [i:] )
		new_domain.extend( d     for d in domain[1:i] )
		new_range .extend( r     for r in range [1:i] )
		
		new_domain.append(Fraction(1))
		new_range .append(image_of_one)
		self._translate_near_zero(new_range)
		return super()._normalise_breakpoints(new_domain, new_range)
	
	def _validate_multiplication(self, other):
		pass
	
	@staticmethod
	def _uncycle(list):
		#Slightly awkward hack here...
		for i in range(len(list) - 1):
			if list[i] > list[i+1]:
				list[i+1] += 1
	
	@classmethod
	def rotation(cls, angle):
		return cls([0, 1], [angle, angle + 1])
	
	def image(self, x):
		"""Take in anything; spew out something in [0, 1)"""
		return super().image(x % 1) % 1
	
	def inverse_image(self, y):
		"""Take in anything; spew out something in [0, 1)"""
		#shift by an integer so that y + n \in [r, r+1]
		r = self.range[0]
		offset = (y - r) % 1
		return super().inverse_image(r + offset)
	
	def restriction_of_range(self, t0, t1):
		d0, d1 = ends(self.range)
		offset = ceil(d0 - t0)
		t0 += offset
		t1 += offset
		domain, range = super().restriction_of_range(t0, t1, raw=True)
		range = [r - offset for r in range]
		return linear_superclass(self)(domain, range)
		
class CPL2(CPLMap, PL2):
	pass