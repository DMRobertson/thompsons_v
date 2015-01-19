from test import setup_script
setup_script(__file__)

import os

from thompson import *
from thompson.examples import *
from thompson.word import format as fmt

"""Uses the `dot` program from the graphviz set of tools to render automorphisms. The `dot` program must be accessible from a folder on the PATH."""

def format(word):
	return fmt(word).replace(' ', '_')

def render_basis(aut, name):
	'''only works for quasinormal'''
	basis = aut.quasinormal_basis
	with open(name + '.dot', 'wt') as f:
		f.write('strict digraph G {\n')
		
		for i in range(1, basis.signature.alphabet_size+1):
			f.write('\tsubgraph cluster{} {{\n'.format(i))
			# f.write('\t\tstyle=filled;\n'.format(i))
			# f.write('\t\tcolor=lightgrey;\n'.format(i))
			f.write('\t\tx{};\n'.format(i))
			
			for j, word in enumerate(basis):
				if word[0] == i:
					f.write('\t\t')
					for k in range(1, len(word)):
						head, _ = word.split(k)
						f.write(format(head) + ' -> ')
					f.write(format(word) + ';\n')
					j2 = basis.index(aut.image(word, inverse=True))
					f.write('\t\t{}[label = "{};{}"];\n'.format(format(word), j+1, j2+1))
			f.write('\t}\n')
		f.write('}\n')
	
	os.system("dot -T png -o {}.png {}.dot".format(name, name))

if __name__ == '__main__':
	render_basis(periodic_QNB_206, 'periodic_QNB_206')
	render_basis(periodic_QNB_344, 'periodic_QNB_344')
	render_basis(cyclic_order_six, 'cyclic_order_six')