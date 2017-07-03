from fractions import Fraction
from .util import all_satisfy, ends, grad, increasing_sequence, int_power_of_two, lerp, pairwise

__all__ = ["PLMap", "PL2"]

class PLMap:
	r"""
	Let :math:`I` and :math:`J` be closed intervals.
	We define :math:`\operatorname{PL}_{S, G}(I, J)` to be the set of functions :math:`\alpha \colon I \to J` with the following properties.
	
	- :math:`\alpha` is an orientation-preserving bijection.
	- :math:`\alpha` is piecewise-linear, with finitely many linear segments.
	- The coordinates of each breakpoint of :math:`\alpha` belong to :math:`S`.
	- The gradients of each linear section of :math:`\alpha` belong to :math:`G`.
	
	We do not require that the endpoints of :math:`I` or :math:`J` belong to :math:`S`.
	
	This class represents functions of this form where :math:`S = \mathbb{Q}` and :math:`G = \mathbb{Q}_{>0}`.
	
	:ivar domain: see :meth:`__init__`.
	:ivar range: see :meth:`__init__`.
	:ivar gradients: The ith entry of this list is the gradient of the ith linear segment.
		
	"""
	def __init__(self, domain, range):
		r"""
		Create a new ``PLMap`` given the breakpoints' coordinates.
		Breakpoint lists are normalised in memory: redudant breakpoints (where the the gradient does not change) are removed.
		Coordinates are provided via two lists *domain* and *range* of :class:`~py3:fractions.Fraction` s.
		
		:raises ValueError: if ``len(domain) != len(range)``.
		:raises ValueError: if ``len(domain) < 2``.
		:raises ValueError: if *domain* or *range* are not increasing sequences.
		:raises ValueError: if *domain* or *range* contain invalid breakpoints. This requirement does **not** apply the first or last element of each list.
		:raises ValueError: if *domain* or *range* describe linear segments with invalid gradients.
		
		.. doctest::
			
			>>> PLMap([0, 1], [0, Fraction(1, 2), 1])
			Traceback (most recent call last):
			...
			ValueError: Domain and range lengths differ
			>>> PLMap([], [])
			Traceback (most recent call last):
			...
			ValueError: Domain must be defined by at least two points
			>>> PLMap([0, 0], [1, 0])
			Traceback (most recent call last):
			...
			ValueError: domain is not an increasing sequence
			>>> PLMap([0, 1], [1, 0])
			Traceback (most recent call last):
			...
			ValueError: range is not an increasing sequence
		"""
		domain = [Fraction(d) for d in domain]
		range  = [Fraction(r) for r in range ]
		self._validate_inputs(domain, range)
		self.domain, self.range, self.gradients = self._normalise_breakpoints(domain, range)

	@classmethod
	def _validate_inputs(cls, domain, range):
		if len(domain) != len(range):
			raise ValueError("Domain and range lengths differ")
		if len(domain) < 2:
			raise ValueError("Domain must be defined by at least two points")
		
		for name, list in (("domain", domain), ("range", range)):
			if not increasing_sequence(list):
				raise ValueError(name + " is not an increasing sequence")
			if not all_satisfy(list[1:-1], cls._validate_breakpoint):
				raise ValueError(name + " contains an invalid breakpoint for " + cls.__name__)
			if not (isinstance(list[ 0], (Fraction)) and
			        isinstance(list[-1], (Fraction))):
				raise ValueError(name + "'s endpoints should be rational")
	
	@classmethod
	def _normalise_breakpoints(cls, domain, range):
		"""
		.. doctest::
		
			>>> x = PLMap(range(5), range(5))
			>>> x.domain
			(Fraction(0, 1), Fraction(4, 1))
			>>> x.range
			(Fraction(0, 1), Fraction(4, 1))
		"""
		initial_gradient = grad(domain[0], domain[1], range[0], range[1])
		if not cls._validate_gradient(initial_gradient):
			raise ValueError("Invalid gradient")
		gradients = [ initial_gradient ]
		
		i = 2
		previous_gradient = initial_gradient
		while i < len(domain):
			next_gradient = grad(domain[i-1], domain[i], range[i-1], range[i])
			if previous_gradient == next_gradient:
				del domain[i-1]
				del range[i-1]
			elif not cls._validate_gradient(initial_gradient):
				raise ValueError("Invalid gradient")
			else:
				gradients.append(next_gradient)
				i += 1
				previous_gradient = next_gradient
		
		return tuple(domain), tuple(range), tuple(gradients)
	
	@classmethod
	def _validate_gradient(cls, gradient):
		return gradient > 0 and isinstance(gradient, Fraction)
	
	@classmethod
	def _validate_breakpoint(cls, breakpoint):
		return isinstance(breakpoint, Fraction)
	
	def __iter__(self):
		r"""Iterating over a ``PLMap`` yields its breakpoints."""
		yield from zip(self.domain, self.range)
		
	@classmethod
	def from_aut(cls, aut):
		"""Creates a new PLMap using a :class:`~thompson.homomorphism.Homomorphism`."""
		domain = []
		range = []
		for segment in aut.pl_segments():
			domain.append(segment['xstart'])
			range.append(segment['ystart'])
		domain.append(segment['xend'])
		range.append(segment['yend'])
		return cls(domain, range)
	
	@classmethod
	def from_stream(cls, stream):
		...
	
	def save_to_file(self, filename):
		...
	
	def __mul__(self, other):
		"""
		Postcompose on the right.
		"""
		if not isinstance(other, type(self)):
			return NotImplemented
		self._validate_multiplication(other)
		domain, range = self._compute_product_breakpoints(other)
		return type(self)(domain, range)
	
	def _validate_multiplication(self, other):
		if ends(self.range) != ends(other.domain):
			raise ValueError("left multiplicand's range is different from right multiplicand's domain")
	
	def _compute_product_breakpoints(self, other):
		domain = set(self.domain)
		for d in other.domain:
			domain.add(self.inverse_image(d))
		domain = sorted(domain)
		range = [other.image(self.image(d)) for d in domain]
		return domain, range
	
	def image(self, x):
		"""Where does the current PLMap send the point :math:`x`?
		
		:raises ValueError: if :math:`x` is not in the current map's domain.
		"""
		if not self.domain[0] <= x <= self.domain[-1]:
			raise ValueError(str(x) + " is not in the domain")
		for (x0, y0), (x1, y1) in pairwise(self):
			if x0 <= x <= x1:
				return lerp(x, x0, x1, y0, y1)
	
	def inverse_image(self, y):
		"""Where is mapped by the current PLMap to the point :math:`y`?
		
		:raises ValueError: if :math:`y` is not in the current map's range.
		"""
		if not self.range[0] <= y <= self.range[-1]:
			raise ValueError(str(y) + " is not in the range")
		for (x0, y0), (x1, y1) in pairwise(self):
			if y0 <= y <= y1:
				return lerp(y, y0, y1, x0, x1)
	
	def __eq__(self, other):
		return self.domain == other.domain and self.range == other.range
	
	def is_identity(self):
		return len(self.domain) == len(self.range) == 2 and domain == range
	
	def __invert__(self):
		return type(self)(self.range, self.domain)
	
	def __xor__(self, other):
		if ends(other.domain) != ends(self.range) or ends(self.domain) != ends(other.range):
			raise ValueError("Mismatching domains and ranges")
		return ~other * self * other
	
	def commutes(self, other):
		return self * other == other * self
	
	def centralise_in_F(*params):
		...

class PL2(PLMap):
	r"""
	:math:`\operatorname{PL}_2` is shorthand for the set of ``PLMap`` s with :math:`S = \mathbb{Z}[1/2]` and :math:`G = 2^{\mathbb Z}`.
	
	.. doctest::
	
		>>> PL2([0, Fraction(1, 3)], [1, Fraction(5, 3)])
		<plmaps.plmap.PL2 object at 0x...>
		>>> PL2([0, Fraction(1, 2), 1], [1, Fraction(5, 3), 2])
		Traceback (most recent call last):
		...
		ValueError: range contains an invalid breakpoint
		>>> PL2([0, 1], [0, 3])
		Traceback (most recent call last):
		...
		ValueError: Invalid gradient
	
	"""
	@classmethod
	def _validate_gradient(cls, gradient):
		if not super()._validate_gradient(gradient):
			return False
		if gradient > 1:
			return gradient.denominator == 1 and int_power_of_two(gradient.numerator  )
		else:
			return gradient.numerator   == 1 and int_power_of_two(gradient.denominator)
	
	@classmethod
	def _validate_breakpoint(cls, breakpoint):
		return super()._validate_breakpoint(breakpoint) and int_power_of_two(breakpoint.denominator)
