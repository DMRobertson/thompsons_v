"""
.. testsetup::
	
	from thompson.examples import *
"""
import os.path
import pkg_resources
import string

from ..generators   import Generators
from ..homomorphism import Homomorphism
from ..automorphism import Automorphism
from ..mixed        import MixedAut
from ..periodic     import PeriodicAut
from ..infinite     import InfiniteAut

from . import random
from .random import *

__all__ = random.__all__

#TODO. Allow powers in the definition of words e.g. a1^4?

def add_module_attribute(name, value):
	globals()[name] = value
	# print('loading', name)
	__all__.append(name)

remove_whitespace = str.maketrans('', '', string.whitespace)

def read_examples():
	#1. Read in examples.
	files = pkg_resources.resource_listdir("thompson", "examples")
	for filename in files:
		name, ext = os.path.splitext(filename)
		if ext != '.aut':
			continue
		path = pkg_resources.resource_filename("thompson.examples", filename)
		aut = Automorphism.from_file(path)
		add_module_attribute(name, aut)
	
	#2. Extract any free factors from mixed periodic/infinite automorphisms.
	free_factors = pkg_resources.resource_filename("thompson.examples", "free_factors.txt")
	with open(free_factors, encoding='utf-8') as f:
		for line in f:
			name = line.lower().strip()
			aut = globals()[name]
			p, i = aut.free_factors()
			add_module_attribute(name + '_p', p)
			add_module_attribute(name + '_i', i)
			p.__doc__ = "A purely periodic free factor, extracted from {}.".format(name)
			i.__doc__ = "A purely infinite free factor, extracted from {}.".format(name)
	
	#3. If any examples have more than one name, deal with that next.
	aliases = pkg_resources.resource_filename("thompson.examples", "aliases.txt")
	with open(aliases, encoding='utf-8') as f:
		for line in f:
			alias, name = line.translate(remove_whitespace).split('=')
			add_module_attribute(alias, globals()[name])

read_examples()