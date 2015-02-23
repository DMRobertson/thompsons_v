"""
.. testsetup::
	from thompson.examples import *
"""

from random import randint, shuffle 
from copy   import copy

from ..word         import Signature, Word
from ..generators   import Generators
from ..automorphism import Automorphism
from ..mixed        import MixedAut
from ..periodic     import PeriodicAut
from ..infinite     import InfiniteAut

__all__ = ['random_signature', 'random_simple_word', 'random_basis',
	'random_automorphism',   'random_periodic_automorphism',   'random_infinite_automorphism',
	'random_automorphism_in_F', 'random_automorphism_in_T',
	'random_conjugate_pair', 'random_conjugate_periodic_pair', 'random_conjugate_infinite_pair',
	'random_powers', 'random_power_bounds', 'random_power_conjugate_pair']

def needs_defaults(undec):
	def decd(signature=None, num_expansions=None):
		if signature is None:
			signature = random_signature()
		if num_expansions is None:
			num_expansions = randint(1, 5)
		return undec(signature, num_expansions)
	decd.__doc__ = undec.__doc__
	return decd

def random_signature():
	r"""Randomly generates a :class:`~thompson.word.Signature` :math:`(n, r)` for use in the functions below. The values of :math:`n` and :math:`r` are chosen (uniformly) at random from :math:`n \in \{2, 3, 4\}` and :math:`r \in \{1, 2, 3, 4, 5\}`, respectively.
	
	.. note:: This function is used to select a random signature when no *signature* argument is provided to the following random functions.
	"""
	return Signature(randint(2, 4), randint(1, 5))

def random_simple_word(signature=None):
	r"""Randomly generates a :meth:`simple <thompson.word.Word.is_simple>` :class:`~thompson.word.Word` belonging to the algebra with the given *signature*. The word consists of an :math:`x_i` followed by 0--15 descendant operators :math:`\alpha_j`. The length and the index of each :math:`\alpha_j` is chosen (uniformly) at random.
	
		>>> random_simple_word().is_simple()
		True
	"""
	signature = signature or random_signature()
	letters = []
	base = randint(1, signature.alphabet_size)
	letters.append(base)
	
	num_alphas = randint(0, 15)
	for i in range(num_alphas):
		letters.append(-randint(1, signature.arity))
	return Word(letters, signature, preprocess=False)

@needs_defaults
def random_basis(signature, num_expansions):
	r"""Randomly generates a basis for the algebra with the given *signature* :math:`(n,r)`. The basis is generated by expanding the :meth:`~thompson.generators.Generators.standard_basis` *num_expansions* times. The expansion point is chosen (uniformly) at random each time. If *num_expansions* is not provided, a value from :math:`\{1, 2, 3, 4, 5\}` is chosen (uniformly) at random.
	
		>>> random_basis().is_basis()
		True
	
	.. note:: This does not generate *bases* uniformly at random. For instance, take :math:`V_{n,r} = V_{2,1}` and let *num_expansions = 3*. The first expansion always gives the basis :math:`[x\alpha_1, x\alpha_2]`. Expanding this twice produces give six bases, one of which appears twice. (To see this, enumerate the rooted binary trees with four leaves.)
	"""
	basis = Generators.standard_basis(signature)
	for _ in range(num_expansions):
		i = randint(0, len(basis) - 1)
		basis.expand(i)
	return basis

@needs_defaults
def random_automorphism(signature, num_expansions):
	"""Randomly generates an automorphism for the algebra with the given signature. Two bases *domain* and *range* are generated by :func:`random_basis`; then *range* is :func:`randomly shuffled <py3:random.shuffle>`. An automorphism is returned which maps the elements of *domain* to those of *range* in order.
	
	.. note:: The bases may be reduced when the automorphism is created (if they contain redundancy).
	"""
	domain = random_basis(signature, num_expansions)
	range  = random_basis(signature, num_expansions)
	shuffle(range)
	return Automorphism(domain, range)

@needs_defaults
def random_periodic_automorphism(signature, num_expansions):
	"""Randomly generates an automorphism for the algebra with the given signature. A basis *domain* is generated by :func:`random_basis`. A copy *range* is made of *domain*, and is then :func:`randomly shuffled <py3:random.shuffle>`. An automorphism is returned which maps the elements of *domain* to those of *range* in order.
	
	.. note:: The bases may be reduced when the automorphism is created (if they contain redundancy.
	"""
	domain = random_basis(signature, num_expansions)
	range  = copy(domain)
	shuffle(range)
	return Automorphism(domain, range)

@needs_defaults
def random_infinite_automorphism(signature, num_expansions):
	"""Randomly generates an infinite automorphism---either by chance, or by extracting an infinite :meth:`~thompson.mixed.MixedAut.free_factor`."""
	if num_expansions == 1:
		#cannot generate an infinite aut with only one expansion
		num_expansions = 2
	phi = None
	while not isinstance(phi, InfiniteAut):
		phi = random_automorphism(signature, num_expansions)
	return phi

@needs_defaults
def random_conjugate_pair(signature, num_expansions):
	r"""Calls :func:`random_automorphism` to create two automorphisms :math:`\psi` and :math:`\rho`. Returns the pair :math:`(\psi, \rho^{-1}\psi\rho)`, which are conjugate by definition."""
	psi = random_automorphism(signature, num_expansions)
	rho = random_automorphism(signature)
	phi = ~rho * psi * rho
	return psi, phi

@needs_defaults
def random_conjugate_periodic_pair(signature, num_expansions):
	r"""The same as :func:`random_conjugate_pair`, except :math:`\psi` is chosen to be a :func:`random_periodic_automorphism`."""
	psi = random_periodic_automorphism(signature, num_expansions)
	rho = random_automorphism(signature)
	phi = ~rho * psi * rho
	return psi, phi

@needs_defaults
def random_conjugate_infinite_pair(signature, num_expansions):
	r"""The same as :func:`random_conjugate_pair`, except :math:`\psi` is chosen to be a :func:`random_infinite_automorphism`."""
	psi = random_infinite_automorphism(signature, num_expansions)
	rho = random_automorphism(signature)
	phi = ~rho * psi * rho
	return psi, phi

def random_powers():
	a = randint(1, 5)
	b = randint(1, 5)
	if randint(0, 1): a *= -1
	if randint(0, 1): b *= -1
	return a, b

#For cheating in doctests. The maximum moduli of a and b that can be produced by the method above
random_power_bounds = (5, 5) 

@needs_defaults
def random_power_conjugate_pair(signature, num_expansions):
	r"""Generates :math:`\psi, \rho` using :func:`random_automorphism` and random powers :math:`a, b`. Returns the tuple :math:`(\Psi, \Phi) = (\psi^b, \rho^{-1}\psi^a\rho)`. Note that :math:`\Psi^a` is conjugate to :math:`\Phi^b` via :math:`\rho`."""
	unfinished = True
	while unfinished:
		psi = random_automorphism(signature, num_expansions)
		a, b = random_powers()
		psi_ab = psi ** (a * b)
		unfinished = psi.order < float('inf') and (a * b) % psi.order == 0
	rho = random_automorphism(signature)
	phi = (~rho) * (psi ** a) * rho
	psi_b = psi **  b
	return psi_b, phi

@needs_defaults
def random_automorphism_in_F(signature, num_expansions):
	"""Generates a random element of :math:`F_{n,r}`."""
	domain = random_basis(signature, num_expansions)
	range = random_basis(signature, num_expansions)
	return Automorphism(domain, range)

@needs_defaults
def random_automorphism_in_T(signature, num_expansions):
	"""Generates a random element of :math:`T_{n,r} \setminus F_{n,r}`."""
	domain = random_basis(signature, num_expansions)
	range = random_basis(signature, num_expansions).cycled(randint(1, len(domain)))
	return Automorphism(domain, range)