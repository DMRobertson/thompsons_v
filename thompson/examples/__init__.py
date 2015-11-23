"""This module provides a number of explicit examples for use in doctests. Additionally, functions to generate random automorphisms are provided.

.. testsetup::
	
	from thompson.examples import *
"""

import os.path
import pkg_resources
import string
from itertools import chain

from ..automorphism import Automorphism
from ..generators   import Generators
from ..word         import Word
from .              import random
from .random        import *

__all__ = random.__all__ + ['available_examples', 'show_examples',
  'load_example', 'load_example_pair',
  'standard_generator'
  ]

#TODO. Allow powers in the definition of words e.g. a1^4?

cache = {}
aliases = None

def load_example(name):
	"""Loads the example with the given *name* from disk. A corresponding :class:`~thompson.automorphism.Automorphism` instance is created and returned. The results are cached, so call this method as often as you like."""
	try:
		return cache[name]
	except KeyError:
		pass
	
	if aliases is None:
		retreive_aliases()
	
	try:
		alt = aliases[name]
	except KeyError:
		path = pkg_resources.resource_filename("thompson.examples", name + '.aut')
		aut = Automorphism.from_file(path)
		cache[name] = aut
	else:
		cache[name] = aut = load_example(alt)
	return aut

def retreive_aliases():
	global aliases
	aliases = dict()
	remove_whitespace = str.maketrans('', '', string.whitespace)
	path = pkg_resources.resource_filename("thompson.examples", "aliases.txt")
	with open(path, encoding='utf-8') as f:
		for line in f:
			alias, name = line.translate(remove_whitespace).split('=')
			aliases[alias] = name

def load_example_pair(name):
	"""Loads a pair of examples, ``*name*_psi`` and ``*name*_phi``.
	
	:rtype: a 2-tuple of automorphisms.
	"""
	return load_example(name + '_psi'), load_example(name + '_phi')

def available_examples():
	"""Returns an iterator yielding the names of the examples that are provided with the package. (Note that :ref:`the full list is provided <examples-table>` in this documentation.)
	
		>>> list(available_examples())[:4]
		['alphabet_size_two', 'arity_four', 'arity_three_order_inf', 'bleak_alpha']
	"""
	files = pkg_resources.resource_listdir("thompson", "examples")
	for filename in sorted(files):
		name, ext = os.path.splitext(filename)
		if ext.lower() == '.aut':
			yield name

def show_examples():
	"""Prints the list of :func:`available_examples`. This is intended for use when interacting with the interpreter."""
	for i, name in enumerate(available_examples()):
		print("{:>3}: {}".format(i+1, name))

def load_all_examples():
	"""Loads (and processes) **all** examples provided by the package. Returns a dictionary whose keys are the example names and whose values are the loaded automorphisms.
	
	.. note::
		To discourage use of this function, it is not available when importing ``*`` from ``thompson.examples``.
		Instead it must be explicitly imported with ``from thompson.examples import load_all_examples``.
	
	.. warning::
		Some of the examples are slow to process, so calling this could take a long time.
	"""
	retreive_aliases()
	for key in chain(available_examples(), aliases):
		load_example(key)
	return cache

def standard_generator(n=0):
	"""Produces the standard generator :math:`X_n` of :math:`F` as described in [CFP96]_."""
	domain  = Generators((2, 1))
	path = Word('x', (2, 1))
	for i in range(n - 1):
		domain.append(path.alpha(1))
		path = path.alpha(2)
	
	#At this stage, domain is equal to the intersection of domain and range
	range_basis = domain.copy()
	
	domain.append(path.alpha(1))
	domain.append(path.extend('a2 a1'))
	domain.append(path.extend('a2 a2'))
	
	range_basis.append(path.extend('a1 a1'))
	range_basis.append(path.extend('a1 a2'))
	range_basis.append(path.alpha(2))
	
	return Automorphism(domain, range_basis)

