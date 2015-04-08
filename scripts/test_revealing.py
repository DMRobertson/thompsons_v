from test import setup_script
setup_script(__file__)

"""Checks all the examples in thompson. Those which are not revealing wrt the minimal tree pair (or alternatively the tree pair corresp to the quasinormal_basis) are listed and rendered."""

from thompson import *
from thompson.examples import load_all_examples
import tempfile, os, os.path, traceback, shutil

# BASES = 'MINIMAL'
BASES = 'QUASINORMAL'

if __name__ == '__main__':
	examples = load_all_examples()
	with tempfile.TemporaryDirectory() as dir:
		os.chdir(dir)

		for name in sorted(examples):
			ex = examples[name]
			if BASES == 'MINIMAL':
				domain = None
			elif BASES == 'QUASINORMAL':
				domain = ex.quasinormal_basis.minimal_expansion_for(ex)
			revealing, data = ex.test_revealing(domain)
			if not revealing:
				print(name, ex.signature, data)
				if input('Render? (y/n) ').strip().lower() == 'y':
					try:
						ex.render(name, name=name.replace('_', r'\_'))
					except Exception:
						traceback.print_exc()
						print('continuing to next example')
		input('Press enter to remove temporary files.\nDrawings currently open will not be removed.')
	print('\nTemporary files removed from', dir)