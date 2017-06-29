from fractions import Fraction
from .util import all_satisfy, grad, increasing_sequence, int_power_of_two, lerp, pairwise

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
	:ivar gradients: see :meth:`__init__`.
		
	"""
	def __init__(self, domain, range):
		r"""
		Create a new ``PLMap`` given the breakpoints' coordinates.
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
		if len(domain) != len(range):
			raise ValueError("Domain and range lengths differ")
		if len(domain) < 2:
			raise ValueError("Domain must be defined by at least two points")
		for name, list in (("domain", domain), ("range", range)):
			if not increasing_sequence(list):
				raise ValueError(name + " is not an increasing sequence")
			if not all_satisfy(list[1:-1], self._validate_breakpoint):
				raise ValueError(name + "contains an invalid breakpoint for " + self.__class__.__name__)
			if not (isinstance(list[ 0], (Fraction, int)) and
			        isinstance(list[-1], (Fraction, int))):
				raise ValueError(name + "'s endpoints should be rational")
		
		#TODO: normalise domain and range: delete redundant breakpoints
		
		self.domain = tuple(Fraction(d) for d in domain)
		self.range  = tuple(Fraction(r) for r in range)
		self.gradients = tuple(
			grad(start[0], end[0], start[1], end[1]) for start, end in pairwise(self)
		)
		if not all_satisfy(self.gradients, self._validate_gradient):
			raise ValueError("Invalid gradient")
	
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
		
	def image(self, x):
		"""Where does the current PLMap send :math:`x`?
		
		:raises ValueError: if :math:`x` is not in the current map's domain.
		"""
		if not self.domain[0] <= x <= self.domain[-1]:
			raise ValueError(str(x) + " is not in the domain")
		for (x0, x1), (y0, y1) in pairwise(self):
			if start[0] <= x <= end[0]:
				return lerp(x, x0, x1, y0, y1)
	
	@classmethod
	def _validate_gradient(cls, gradient):
		return gradient > 0 and isinstance(gradient, Fraction)
	
	@classmethod
	def _validate_breakpoint(cls, breakpoint):
		return isinstance(breakpoint, Fraction)
	
	def __iter__(self):
		r"""Iterating over a ``PLMap`` yields its breakpoints."""
		yield from zip(self.domain, self.range)

class PL2(PLMap):
	r"""
	:math:`\operatorname{PL_2}` is shorthand for the set of ``PLMap`` s with :math:`S = \mathbb{Z}[1/2]` and :math:`G = 2^{\mathbb Z}`.
	
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
			return gradient.denominator == 1 and int_power_of_two(gradient.numerator)
		else:
			return gradient.numerator == 1 and int_power_of_two(gradient.denominator)
	
	@classmethod
	def _validate_breakpoint(cls, breakpoint):
		return super()._validate_breakpoint(breakpoint) and int_power_of_two(breakpoint.denominator)
