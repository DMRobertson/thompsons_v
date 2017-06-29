from fractions import Fraction
from .util import all_satisfy, grad, pairwise

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
		"""
		if len(domain) != len(range):
			raise ValueError("Domain and range lengths differ")
		if len(domain) < 2:
			raise ValueError("Domain must be defined by at least two points")
		for list, name in (("domain", domain), ("range", range)):
			if not increasing_sequence(list):
				raise ValueError(name + " is not an increasing sequence")
			if not all_satisfy(list[1:-1], self._validate_breakpoint):
				raise ValueError(name + "contains an invalid breakpoint for " + self.__class__.__name__)
		self.domain = tuple(domain)
		self.range  = tuple(range)
		self.gradients = tuple(
			grad(start, end) for start, end in pairwise(self)
		)
		if not all_satisfy(self.gradients, self._validate_gradient):
			raise ValueError("Invalid gradient")
	
	@staticmethod
	def _validate_gradient(gradient):
		return gradient > 0 and isinstance(gradient, Fraction)
	
	@staticmethod
	def _validate_breakpoint(breakpoint):
		return isinstance(breakpoint, Fraction)
	
	def __iter__(self):
		r"""Iterating over a ``PLMap`` yields its breakpoints."""
		yield from zip(self.domain, self.range)

class PL2(PLMap):
	r"""
	:math:`\operatorname{PL_2}` is shorthand for the set of ``PLMap`` s with :math:`S = \mathbb{Z}[1/2]` and :math:`G = 2^{\mathbb Z}`.
	"""
	@staticmethod
	def _validate_gradient(gradient):
		if not super()._validate_gradient(gradient):
			return False
		if gradient < 1:
			gradient = 1 / gradient
		return gradient.denominator == 1 and int_power_of_two(gradient.numerator)
	
	@staticmethod
	def _validate_breakpoint(breakpoint):
		return super()._validate_breakpoint(breakpoint) and int_power_of_two(gradient.denominator)
