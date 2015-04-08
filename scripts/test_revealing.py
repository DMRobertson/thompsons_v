from test import setup_script
setup_script(__file__)

"""A place to play around with the package."""

from thompson import *
from thompson.examples import load_all_examples
import tempfile, os, os.path, traceback, shutil

if __name__ == '__main__':
	examples = load_all_examples()
	with tempfile.TemporaryDirectory() as dir:
		os.chdir(dir)

		for name in sorted(examples):
			ex = examples[name]
			revealing, data = ex.test_revealing()
			if not revealing:
				print(name, ex.signature, data)
				if input('Render? (y/n) ').strip().lower() == 'y':
					try:
						ex.render(name, name=name.replace('_', r'\_'))
					except Exception:
						traceback.print_exc()
						print('continuing to next example')
	print('\nTemporary files were stored in', dir, 'and some may still be present.')