from test import setup_script
from enum import Enum
setup_script(__file__)

"""Checks all the examples in thompson. Those which are not revealing wrt the minimal tree pair (or alternatively the tree pair corresp to the quasinormal_basis) are listed and rendered."""

from thompson import *
from thompson.examples import load_all_examples
import tempfile, os, os.path, traceback, shutil, re

BASES = 'MINIMAL'
# BASES = 'QUASINORMAL'

class BasisOptions(Enum):
	minimal = 1
	quasinormal = 2

desc = {
	BasisOptions.minimal: 'Use the minimal tree pair describing the automorphism.',	
	BasisOptions.quasinormal: 'Use the tree pair corresponding to the minimal expansion of the quasinormal basis.'
}

def choose_from_enum(enum, desc = None):
	if desc is None:
		desc = {}
	num_choices = len(enum)
	print('Choices are:')
	for i, x in enumerate(enum):
		print('\t[{}] {}. {}'.format(i + 1, x.name, desc.get(x, '')))

		ans = None
	while ans is None:
		try:
			ans = int(input('Please make a choice by entering an integer in [1--{}]: '.format(
			  num_choices
			)))
		except Exception as e:
			print(e)
		else:
			if not (1 <= ans <= num_choices):
				ans = None
	return enum(ans)

if __name__ == '__main__':
	print('This script looks through all the examples of thompson for automorphisms with non-revealing tree pairs.')

	basis_option = choose_from_enum(BasisOptions, desc)
	examples = load_all_examples()
	with tempfile.TemporaryDirectory() as dir:
		os.chdir(dir)
		print('\nCreated temporary directory:', dir)
		
		for name in sorted(examples):
			ex = examples[name]
			if basis_option == BasisOptions.minimal:
				domain = None
			elif basis_option == BasisOptions.quasinormal:
				domain = ex.quasinormal_basis.minimal_expansion_for(ex)
			revealing, data = ex.test_revealing(domain)
			if not revealing:
				print('The automorphism {} of V_{} is not revealing (consider {}).'.format(
				  name, ex.signature, data))
				if input('Render? (y/n) ').strip().lower() == 'y':
					try:
						ex.render(name, name=name.replace('_', r'\_'))
					except Exception:
						traceback.print_exc()
						print('Continuing to next example.')
		input('Press enter to remove temporary files.\nDrawings currently open will not be removed.')
	print('\nTemporary files removed from', dir)
