"""
.. testsetup::
	from thompson.automorphism import *
"""

class Automorphism:
	r"""Represents an automorphism of :math:`V_{n,r}` by specifying two bases. This class keeps track of the mapping between bases.
	
	:ivar arity: :math:`n`, the number of operators :math:`\alpha_i`.
	:ivar alphabet_size: :math:`r`, the number of letters :math:`x_i`.
	:ivar domain: a :class:`generating set <Generators>` of preimages 
	:ivar range: a :class:`generating set <Generators>` of images.
	"""
	
	def __init__(self, arity, alphabet_size, domain, range):
		"""Creates an automorphism, given the parameters :math:`n` and :math:`r`. Two bases *domain* and *range* are given. The automorphism maps elements as follows:
		
			.. math:: \text{domain}_i \mapsto \text{range}_i \mapsto \text{range}_{\text{perm}(i)}
		
		:raises ValueError: if the bases are of different sizes.
		:raises ValueError: if the bases have different arities or alphabet sizes.
		:raises ValueError: if either basis isn't actually a basis, i.e. is not a :meth:`free generating set <Generator.is_free>` or does not :meth:`contract to the standard basis <is_basis>`.
		"""
		
		#The boring checks
		if len(domain) != len(range):
			raise ValueError("Domain basis has {} elements, but range basis has {} elements.".format(
			  len(domain), len(range)))
		
		if not(arity == domain.arity == range.arity):
			raise ValueError("Arities do not match. Expected: {}, Domain: {}, Range: {}".format(
			  arity, domain.arity, range.arity))
		
		if not(alphabet_size == domain.alphabet_size == range.alphabet_size):
			raise ValueError("Alphabet sizes do not match. Expected: {}, Domain: {}, Range: {}".format(
			  alphabet_size, domain.alphabet_size, range.alphabet_size))
		
		#Check to see that the two generating sets given are free
		i, j = domain.test_free()
		if not(i == j == -1):
			raise ValueError("Domain is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		i, j = range.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at indices {} and {}.".format(
			  i, j))
		
		#Check to see that the two generating sets generate all of V_{n,r}
		missing = domain.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Domain does not generate V_{},{}. Missing elements are {}.".format(
			  arity, alphabet_size, [format(x) for x in missing]))
		
		missing = range.test_generates_algebra()
		if len(missing) > 0:
			raise ValueError("Range does not generate V_{},{}. Missing elements are {}.".format(
			  arity, alphabet_size, [format(x) for x in missing]))
		
		self.arity = arity
		self.alphabet_size = alphabet_size
		self.domain = domain
		self.range = range
		
	def to_quasinormal_form(self):
		pass
	
