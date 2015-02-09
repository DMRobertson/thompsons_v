import os, sys
os.chdir('..')
sys.path.append(os.getcwd())

from thompson import Automorphism
import thompson.examples

automorphisms = []

for name, attr in thompson.examples.__dict__.items():
	if not isinstance(attr, Automorphism):
		continue
	doc = attr.__doc__.split('\n', maxsplit=1)[0]
	automorphisms.append((name, doc))

automorphisms.sort()

widths = [max(len(item[i]) for item in automorphisms) for i in range(2)]
print(widths)
format_string = "{}\n\t{}"


os.chdir('docs')

with open('examples_table.txt', 'wt') as f:
	for name, docstring in automorphisms:
		print(format_string.format(name, docstring), file=f)