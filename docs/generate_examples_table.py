import os, sys
os.chdir('..')
sys.path.append(os.getcwd())

from thompson import Automorphism
import thompson.examples

automorphisms = []

for name, aut in thompson.examples.load_all_examples().items():
	doc = aut.__doc__.split('\n\n', maxsplit=1)[0] #stuff is ignored after the first appearance of \n\n
	automorphisms.append((name, doc.strip()))

automorphisms.sort()

widths = [max(len(item[i]) for item in automorphisms) for i in range(2)]
print(widths)
format_string = "{}\n\t{}"


os.chdir('docs')

with open('examples_table.txt', 'wt') as f:
	for name, docstring in automorphisms:
		print(format_string.format(name, docstring), file=f)