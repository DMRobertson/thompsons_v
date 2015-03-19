from test import setup_script
setup_script(__file__)

"""A place to play around with the package."""

from thompson import Automorphism, Word
olga_f = Automorphism.from_file('../thompson/examples/olga_f.aut')

head = r'''\begin{tikzpicture}[level distance=8mm,
	level 1/.style={sibling distance=25mm},
	level 2/.style={sibling distance=1cm},
	level 3/.style={sibling distance=5mm}]
\usetikzlibrary{positioning}
'''

tail = r'''\end{tikzpicture}'''

def dump_tikz_code(aut, filename, domain=None, name=''):
	if domain is None:
		domain = aut.domain
	range = aut.image_of_set(domain)
	with open(filename, 'wt') as f:
		f.write(head)
		basis_to_tree(domain, range, f, name='domain')
		
		y = max(len(w) for w in domain) - 1
		z = max(len(w) for w in domain) - 1
		depth = -min(y, z)/2
		
		coordend = r'\textwidth, ' + str(depth*8) + 'mm)'
		print(r'\draw[->, thick] (0.25' + coordend, '-- node[auto]{', name, r'} (0.35' + coordend, ';\n', file=f)
		basis_to_tree(range, domain, f, name='range', extra=r'xshift=0.6\textwidth')
		f.write(tail)



def basis_to_tree(basis, other, f, name, extra=''):
	f.write(r'\begin{scope}[' + extra + ']\n')
	f.write(r'\coordinate (' + name + ') {} \n')
	depth = 1
	dashed = []
	
	generator = iter(basis.preorder_traversal())
	next(generator)
	for path, leafnum in generator:
		while depth >= len(path):
			depth -= 1
			f.write('\t' * depth)
			if dashed.pop():
				f.write(' edge from parent[dashed] ')
			print('}', file=f)
		
		dashed.append(path not in other and other.is_above(path))
		
		f.write('\t' * (depth) + 'child {')
		if leafnum is not None:
			f.write('node{' + str(leafnum) + '}')
			if dashed.pop():
				f.write(' edge from parent[dashed] ')
			f.write('}')
		else:
			depth += 1
		print(file=f)
	while depth > 1:
		depth -= 1
		f.write('\t' * depth)
		if dashed.pop():
			f.write(' edge from parent[dashed] ')
		print('}', file=f)
	print(';', file=f)
	print(r'\end{scope}', file=f)

dump_tikz_code(olga_f, 'test.tikz')