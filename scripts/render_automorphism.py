from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.word import format

"""A hack which turns a basis into a string of 1s and 0s for use with older code."""

def split_basis(basis):
	assert basis.is_basis()
	strings = []
	for i in range(1, basis.signature.alphabet_size + 1):
		part = Generators(basis.signature, (w for w in basis if w[0] == i))
		print(part)
		strings.append(to_binary_string(part, basis.signature.arity, i))
	return strings

def to_binary_string(basis, arity, root_number):
	if len(basis) == 0:
		return '0'
	
	string = []
	path = (root_number,)
	
	print(format(path))
	for word in basis:
		print('looking for', word)
		# 1. Descend to word
		while len(word) > len(path):
			path = path + (-1,)
			string.append('1')
			print(format(path))
		# 2. Are we at the word we want to be at?
		while word != path:
			path = next(path, arity)
		print(format(path))
		string.append('0')
		path = next(path, arity)
		print(format(path))
		
	return ''.join(string)

def next(path, arity):
	while -path[-1] == arity:
		path = path[:-1]
	return path[:-1] + (path[-1]-1,)

if __name__ == '__main__':
	aut = cyclic_order_six
	print(aut)
	print(split_basis(aut.domain))
	aut = Automorphism.from_file('QNB_size_206.aut')
	print(split_basis(aut.quasinormal_basis()))