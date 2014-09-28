
class Automorphism:
	
	def __init__(self, domain_elmts, range_elmts, perm=None):
		r"""Represents an automorphism of :math:`V_{n,r}` by specifying
		
		- two bases *domain_elmts* and *range_elmts*, and
		- a permutation *perm* of *range_elmts*.
		
		The automorphism maps elements as follows:
		
			.. math:: \text{domain}_i \mapsto \text{range}_i \mapsto \text{range}_{\text{perm}(i)}
		
		:raises ValueError: if the bases are of different sizes
		:raises ValueError: if the bases aren't actually bases (i.e. they are not a free generating set)
		"""
		
		if len(domain_elmts) != len(range_elmts):
			raise ValueError("Domain basis has {} elements, but range basis has {} elements."
			    .format(len(domain_elmts), len(range_elmts)))
		
		#Check to see that the two generating sets given are free
		i, j = domain_elmts.test_free()
		if not(i == j == -1):
			raise ValueError("Domain is not a free generating set. Check elements at positions {} and {}.".format(i, j))
		
		i, j = range_elmts.test_free()
		if not(i == j == -1):
			raise ValueError("Range is not a free generating set. Check elements at positions {} and {}.".format(i, j))
		
		#Check to see that the two generating sets given generate all of V_{n,r}
		if not domain_elmts.generates_algebra(airty, alphabet_size):
			raise ValueError("Domain does not generate V_.")

def quasinormal_basis(phi):
	"""An implementation of Lemma 4.24.1. Given an automorphism :math:`\phi \in G_{n, r}`, we can construct the quasinormal basis for :math:`\phi` by expanding the :meth:`standard basis <~thompson.generators.standard_basis>`.
	"""
	pass

