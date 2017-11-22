from fractions import Fraction
from .util import all_satisfy, ends, fixed_point, grad, gradient_roots_dyadic, increasing_sequence, int_power_of_two, lerp, pairwise

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
	__slots__ = ('domain', 'range', 'gradients')
	
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
			elif not cls._validate_gradient(next_gradient):
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
	
	@classmethod
	def identity(cls, t0, t1):
		return cls([t0, t1], [t0, t1])
	
	def __iter__(self):
		r"""Iterating over a ``PLMap`` yields its breakpoints."""
		yield from zip(self.domain, self.range)
	
	def __len__(self):
		return len(self.domain)
	
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
		check = self._validate_multiplication(other)
		if check is NotImplemented:
			return NotImplemented
		domain, range = self._compute_product_breakpoints(other)
		return type(self)(domain, range)
	
	def _validate_multiplication(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		s0, s1 = ends(self.range)
		d0, d1 = ends(other.domain)
		if not d0 <= s0 < s1 <= d1:
			raise ValueError("left multiplicand's range is not a subset of right multiplicand's domain")
	
	def _compute_product_breakpoints(self, other):
		from .cplmap import CPLMap
		domain = set(self.domain)
		s0, s1 = ends(self.range)
		for d in other.domain:
			if isinstance(self, CPLMap) or s0 < d < s1:
				domain.add(self.inverse_image(d))
		domain = sorted(domain)
		range = [other.image(self.image(d)) for d in domain]
		
		if isinstance(other, CPLMap) or isinstance(self, CPLMap):
			CPLMap._uncycle(range)
	
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
	
	def rgradient_at(self, x):
		for i, (d0, d1) in enumerate(pairwise(self.domain)):
			if d0 <= x < d1:
				return self.gradients[i]
		raise ValueError(str(x) + " is not in the domain")
	
	def __eq__(self, other):
		return self.domain == other.domain and self.range == other.range
	
	def is_identity(self):
		return len(self.domain) == len(self.range) == 2 and self.domain == self.range
	
	def __invert__(self):
		return type(self)(self.range, self.domain)
	
	def __xor__(self, other):
		if ends(other.domain) != ends(self.range) or ends(other.domain) != ends(self.domain):
			raise ValueError("Mismatching domains and ranges")
		return ~other * self * other
	
	def __str__(self):
		domain = ""
		range  = ""
		for d, r in self:
			d = str(d)
			r = str(r)
			n = max(len(d), len(r))
			domain += d.ljust(n) + "  "
			range  += r.ljust(n) + "  "
		return domain[:-2] + "\n" + range[:-2]
	
	def __repr__(self):
		return "<{}: [{}, {}] -> [{}, {}]>".format(
			type(self).__name__, str(self.domain[0]), str(self.domain[-1]), str(self.range[0]), str(self.range[-1]))
	
	def dump(self, short=True):
		cls = "F" if short else "Fraction"
		dterms = [self._format_fraction(d, cls) for d in self.domain]
		rterms = [self._format_fraction(r, cls) for r in self.range ]
		for i in range(len(self)):
			m = max(len(dterms[i]), len(rterms[i]))
			dterms[i] = dterms[i].ljust(m)
			rterms[i] = rterms[i].ljust(m)
		
		output = "{}(\n\t[{}],\n\t[{}]\n)".format(
			type(self).__name__,
			", ".join(dterms),
			", ".join(rterms)
		)
		return output
	
	@staticmethod
	def _format_fraction(fraction, cls="Fraction"):
		if fraction.denominator == 1:
			return str(fraction.numerator)
		else:
			return  cls + "({},{})".format(fraction.numerator, fraction.denominator)
	
	def tikz_path(self):
		return " -- ".join(
			"({}, {})".format(str(d), str(r))
			for d, r in self
		) + ";"
	
	def commutes(self, other):
		return self * other == other * self
	
	def centralise_in_F(*params):
		...
	
	def restriction(self, t0, t1):
		r"""
		Produce a copy of the current PLMap restricted to the interval :math:`[t_0, t_1]`.
		
		:param iterable target: A sequence of at least two integers. The first and last entries are the start and end of the interval onto which we restrict.
		
		:raises ValueError: if ``t0 >= t1``.
		:raises ValueError: if *t0* and *t1* do not describe a subinterval of ``self.domain``.
		"""
		if not t0 < t1:
			raise ValueError("Start point is not less than end point")
		d0, d1 = ends(self.domain)
		# Problems here with circular?
		if not d0 <= t0 < t1 <= d1:
			raise ValueError("Target interval is not within domain")
		
		domain = [t0]
		range  = [self.image(t0)]
		
		for d, r in self:
			if t0 < d < t1:
				domain.append(d)
				range .append(r)
		
		domain.append(t1)
		range .append(self.image(t1))
		from . import CPLMap
		if isinstance(self, CPLMap):
			self._uncycle(range)
		return linear_superclass(self)(domain, range)
	
	def restriction_of_range(self, t0, t1, raw=False):
		"""
		Produce a copy of the current PLMap restricted to the interval :math:`[t_0, t_1]`.
		
		:param iterable target: A sequence of at least two integers. The first and last entries are the start and end of the interval onto which we restrict.
		
		:raises ValueError: if ``t0 >= t1``.
		:raises ValueError: if *t0* and *t1* do not describe a subinterval of ``self.domain``.
		"""
		if not t0 < t1:
			raise ValueError("Start point is not less than end point")
		d0, d1 = ends(self.range)
		# Problems here with circular?
		if not d0 <= t0 < t1 <= d1:
			raise ValueError("Target interval is not within range")
		
		domain = [self.inverse_image(t0)]
		range  = [t0]
		
		for d, r in self:
			if t0 < r < t1:
				domain.append(d)
				range .append(r)
		
		domain.append(self.inverse_image(t1))
		range .append(t1)
		from . import CPLMap
		if isinstance(self, CPLMap):
			self._uncycle(domain)
		if raw:
			return domain, range
		else:
			return linear_superclass(self)(domain, range)
	
	def is_permutation(self):
		return ends(self.domain) == ends(self.range)
	
	def fixed_points(self, raw=False):
		output = []
		for (d0, r0), (d1, r1) in pairwise(self):
			if d0 == r0 and d1 == r1:
				output.append((d0, d1))
			else:
				x = fixed_point(d0, d1, r0, r1)
				if x is not None:
					output.append(x)
		
		if not raw:
			self._normalise_points_and_intervals(output)
		return tuple(output)
	
	@staticmethod
	def _normalise_points_and_intervals(output):
		i = 0
		while i < len(output) - 1:
			if (
				isinstance(output[i  ], Fraction) and
				isinstance(output[i+1], tuple   ) and
				output[i] == output[i+1][0]
			):
				del output[i]
			elif (
				isinstance(output[i  ], tuple) and
				isinstance(output[i+1], Fraction   ) and
				output[i][1] == output[i+1]
			):
				del output[i+1]
				i += 1
			else:
				i += 1
	
	def is_one_bump(self):
		return ends(self.domain) == ends(self.range) == self.fixed_points(raw=True)
	
	def one_bump_test_conjugate_with(self, other, initial_gradient):
		if self.gradients[0] != other.gradients[0] or self.gradients[-1] != other.gradients[-1]:
			return None
		if not (self.is_one_bump() and other.is_one_bump() and ends(self.domain) == ends(other.domain)):
			raise ValueError("Functions are not one-bump, or have different domains")
		
		start = self.domain[0]
		end   = self.domain[-1]
		conj_linear_upto, sources_linear_from  = self.one_bump_linearity_boxes(other, initial_gradient)
		
		# assert conj_linear_upto and sources_linear_from are valid breakpoints
		
		domain = (start, conj_linear_upto)
		range  = (start, start + (conj_linear_upto - start) * initial_gradient)
		candidate = type(self)(domain, range)
		
		while not(
			candidate.domain[-1] >= sources_linear_from and
			candidate.range [-1] >= sources_linear_from
		):
			LHS = self.restriction_of_range(start, candidate.domain[-1])
			candidate = LHS * candidate * ~other
		
		domain = candidate.domain + (end,)
		range  = candidate.range  + (end,)
		try:
			candidate = type(self)(domain, range)
		except ValueError as e:
			#Should raise own exception class instead of checking message
			if ("Invalid gradient" in e.args or
			    "domain contains an invalid breakpoint for PL2" in e.args or 
			    "range contains an invalid breakpoint for PL2" in e.args):
				return None
		
		if self ^ candidate == other:
			return candidate
	
	def one_bump_linearity_boxes(self, other, initial_gradient):
		sources_linear_upto = min(self.domain[ 1], other.domain[ 1])
		sources_linear_from = max(self.domain[-2], other.domain[-2])
		
		if initial_gradient < 1:
			conj_linear_upto = sources_linear_upto
		else:
			conj_linear_upto = start + (sources_linear_upto - start) / initial_gradient
		
		return conj_linear_upto, sources_linear_from
		

class PL2(PLMap):
	r"""
	:math:`\operatorname{PL}_2` is shorthand for the set of ``PLMap`` s with :math:`S = \mathbb{Z}[1/2]` and :math:`G = 2^{\mathbb Z}`.
	As `noted above <plmaps.plmaps.PLMap>`, we don't insist that the endpoints of the domain and range lie in S: we're working with what my thesis calls :math:`\operatorname{PL}_2^\text{rest}` rather than :math:`\operatorname{PL}_2^\text{flat}`.
	
	.. doctest::
	
		>>> PL2([0, Fraction(1, 3)], [1, Fraction(5, 3)])
		<PL2: [0, 1/3] -> [1, 5/3]>
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
	
	def one_bump_test_conjugate(self, other):
		if not (self.is_one_bump() and other.is_one_bump()):
			raise ValueError("Given functions are not one-bump")
		#Must have self.gradients[0] == other.gradients[0]
		#If there exists a conjugator, there must exist a (possibly different) conjugator with initial gradient in [1, A) or (A, 1] (whichever is nonempty)
		bound = self.gradients[0]
		if bound > 1:
			limits = (1, bound)
		else:
			limits = (bound, 1)
		power = limits[0]
		while power <= limits[1]:
			result = self.one_bump_test_conjugate_with(other, power)
			if result is not None:
				return result
			power *= 2
	
	def one_bump_cent_gen(self):
		"""If the current PLMap is a one-bump function :math:`D \to D`, produce an element which generates its centraliser in :math:`\operatorname{PL}(D)`.
		The generator's initial gradient will be above 1 if and only if the the current PLMap's initial gradient is above 1.
		
		:raises ValueError: if the current PLMap is not a one-bump function.
		"""
		if not self.is_one_bump():
			raise ValueError("Given function is not one-bump")
		initial_gradient = self.gradients[0]
		for target_gradient in gradient_roots_dyadic(initial_gradient):
			result = self.one_bump_test_conjugate_with(self, target_gradient)
			if result is not None:
				return result
		return self

def linear_superclass(self):
	for cls in self.__class__.__bases__:
		if cls in {PLMap, PL2}:
			return cls

def glue(*maps, cls=None):
	"""Assume maps is a list of PLMaps, whose domains and ranges can be glued together properly.
	"""
	domain = maps[0].domain + tuple( x
		for map in maps[1:] for x in map.domain[1:])
	range  = maps[0].range  + tuple( x
		for map in maps[1:] for x in map.range [1:])
	if cls is None:
		if all (isinstance(map, PL2) for map in maps):
			cls = PL2
		else:
			cls = PLMap
	
	return cls(domain, range)
