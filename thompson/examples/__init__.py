"""A list of examples for use in doctests."""

import pkg_resources
import re
import string

from ..generators import Generators
from ..automorphism import Automorphism

__all__ = ["aut_from_file"]

#TODO. Have a method that reads homomorphisms in from a file. Writing the Python syntax around all this is boring.
#TODO. Allow powers in the definition of words e.g. a1^4?

extract_signatures = re.compile(r"""
	\( \s* (\d+) [\s,]+ (\d+) \s* \)
	[\s\->]+
	\( \s* (\d+) [\s,]+ (\d+) \s* \)""", re.VERBOSE)

def aut_from_file(filename):
	"""Reads in a file specifying an automorphism and returns a pair *(aut, name)*. Here is an example of the format::
	
	example_4_5
	5
	(2,1)		->	(2,1)
	x1 a1 a1 a1	->	x1 a1 a1
	x1 a1 a1 a2	->	x1 a1 a2 a1
	x1 a1 a2	->	x1 a1 a2 a2
	x1 a2 a1	->	x1 a2 a2
	x1 a2 a2	->	x1 a2 a1
	False
	False
	
	The specific details:
	
	- name
	- number of generators
	- signatures of domain and range
	- list of rules domain -> range
	- Boolean: should we extract a periodic factor?
	- Boolean: should we extract an infinite factor?
	"""
	with open(filename, encoding='utf-8') as f:
		name = f.readline().strip()
		num_generators = int(f.readline())
		
		params = extract_signatures.match(f.readline().strip()).groups()
		d = Generators([int(params[0]), int(params[1])])
		r = Generators([int(params[2]), int(params[3])])
		
		for _ in range(num_generators):
			d_word, r_word = (word.strip() for word in f.readline().split('->'))
			d.append(d_word)
			r.append(r_word)
		extract_periodic = f.readline().strip() == 'True'
		extract_infinite = f.readline().strip() == 'True' 
	
	aut = Automorphism(d, r)
	return aut, name, extract_periodic, extract_infinite

def add_module_attribute(name, value):
	globals()[name] = value
	__all__.append(name)

remove_whitespace = str.maketrans('', '', string.whitespace)

def read_examples():
	files = pkg_resources.resource_listdir("thompson", "examples")
	files = (filename for filename in files if filename.endswith('.aut'))
	for filename in files:
		path = pkg_resources.resource_filename("thompson.examples", filename)
		aut, name, extract_periodic, extract_infinite = aut_from_file(path)
		add_module_attribute(name, aut)
		
		qnf_basis = aut.quasinormal_basis()
		p, i = aut._partition_basis(qnf_basis)
		if extract_periodic:
			add_module_attribute(name + '_p', aut.free_factor(p, infinite=False))
		if extract_infinite:
			add_module_attribute(name + '_i', aut.free_factor(i, infinite=True))
	
	aliases = pkg_resources.resource_filename("thompson.examples", "aliases.txt")
	with open(aliases, encoding='utf-8') as f:
		for line in f:
			alias, name = line.translate(remove_whitespace).split('=')
			add_module_attribute(alias, globals()[name])

read_examples()
