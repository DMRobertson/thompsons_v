"""
.. testsetup::

	from thompson import F, T, load_example, random_automorphism

As an alternative to using the :meth:`~thompson.automorphism.Automorphism.preserves_order` and :meth:`~thompson.automorphism.Automorphism.cycles_order` methods, we provide two objects ``F``, ``T`` which  know how to respond to the ``in`` operator.

	>>> f = random_automorphism(group='F')
	>>> f in F and f in T
	True
	>>> t = load_example('scott_free_alpha')
	>>> t not in F and t in T
	True
"""

from .automorphism import Automorphism

__all__ = ["F", "T"]

class FClass:
	def __str__(self):
		return 'F'
	
	def __contains__(self, key):
		return isinstance(key, Automorphism) and key.preserves_order()
F = FClass()

class TClass:
	def __str__(self):
		return 'T'
	
	def __contains__(self, key):
		return isinstance(key, Automorphism) and key.cycles_order()
T = TClass()

class VClass:
	def __str__(self):
		return 'V'
	
	def __contains__(self, key):
		return isinstance(key, Automorphism)
V = VClass()