"""
.. testsetup::
	
	from thompson.examples import *
"""

import pkg_resources
import string

from ..generators import Generators
from ..automorphism import Automorphism

from . import random
from .random import *

__all__ = random.__all__

#TODO. Allow powers in the definition of words e.g. a1^4?

def add_module_attribute(name, value):
	globals()[name] = value
	__all__.append(name)

remove_whitespace = str.maketrans('', '', string.whitespace)

def read_examples():
	#1. Read in examples.
	files = pkg_resources.resource_listdir("thompson", "examples")
	files = [filename for filename in files if filename.endswith('.aut')]
	for filename in files:
		path = pkg_resources.resource_filename("thompson.examples", filename)
		aut = Automorphism.from_file(path)
		add_module_attribute(filename[:-4], aut)
	
	#2. Extract any free factors which are useful.
	free_factors = pkg_resources.resource_filename("thompson.examples", "free_factors.txt")
	with open(free_factors, encoding='utf-8') as f:
		for line in f:
			details = line.lower().split()
			name = details[0]
			
			aut = globals()[name]
			basis = aut.quasinormal_basis()
			p, i = aut._partition_basis(basis)
			if 'p' in details[1:]:
				add_module_attribute(name + '_p', aut.free_factor(p, infinite=False))
			if 'i' in details[1:]:
				add_module_attribute(name + '_i', aut.free_factor(i, infinite=True))
	
	#3. If any examples have more than one name, deal with that next.
	aliases = pkg_resources.resource_filename("thompson.examples", "aliases.txt")
	with open(aliases, encoding='utf-8') as f:
		for line in f:
			alias, name = line.translate(remove_whitespace).split('=')
			add_module_attribute(alias, globals()[name])

read_examples()