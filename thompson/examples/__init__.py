"""
.. testsetup::
	
	from thompson.examples import *
"""
import os.path
import pkg_resources
import string

from ..automorphism import Automorphism

from . import random
from .random import *

__all__ = random.__all__ + ['load_example', 'load_example_pair']

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

def load_all_examples():
	"""Loads (and processes) **all** examples provided by the package. Returns a dictionary whose keys are the example names and whose values are the loaded automorphisms.
	
	.. note::
		To discourage use of this function, it is not available when importing ``*`` from ``thompson.examples``.
		Instead it must be explicitly imported with ``from thompson.examples import load_all_examples``.
	
	.. warning::
		Some of the examples are slow to process, so calling this could take a long time.
	"""
	files = pkg_resources.resource_listdir("thompson", "examples")
	for filename in files:
		name, ext = os.path.splitext(filename)
		if ext != '.aut':
			continue
		load_example(filename[:-4])
	retreive_aliases()
	for alias in aliases:
		load_example(alias)	
	return cache