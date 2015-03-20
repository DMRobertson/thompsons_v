from test import setup_script
setup_script(__file__)

"""A place to play around with the package."""

from thompson import *

head = r'''\tikzset{
	level distance=8mm,
	level 1/.style={sibling distance=25mm},
	level 2/.style={sibling distance=15mm},
	level 3/.style={sibling distance=10mm},
	level 4/.style={sibling distance=5mm},
	characteristic/.style={red, draw, solid}
}
'''

tail = '' #r'''\end{tikzpicture}'''

def dump_tikz_code(aut, filename, domain=None, name=''):
	"""Saves instructions in the given *filename* for TikZ to draw the automorphism *aut*. The leaves of the domain tree are the elements of *domain*."""
	#1. Setup
	if domain is None:
		domain = aut.domain
	range = aut.image_of_set(domain)
	#Compute X = L(domain \intersect range)
	X = domain.copy()
	X += range
	finished = False
	while not finished:
		i, j = X.test_free()
		if i == j == -1:
			finished = True
		else:
			del X[j]
	X.sort()
	
	with open(filename, 'wt') as f:
		f.write(head)
		basis_to_tree(aut, domain, range, X, f, name='domain')
		
		y = max(len(w) for w in domain) - 1
		z = max(len(w) for w in domain) - 1
		depth = -min(y, z)/2
		
		coordend = r'\textwidth, ' + str(depth*8) + 'mm)'
		print(r'\draw[->, thick] (0.25' + coordend, '-- node[auto]{$', name, r'$} (0.35' + coordend, ';\n', file=f)
		
		basis_to_tree(aut, range, domain, X, f, name='range', extra=r'xshift=0.6\textwidth')
		f.write(tail)

def basis_to_tree(aut, basis, other, intersection, f, name, extra=''):
	print(r'\begin{scope}[', extra, ']', file=f)
	print(r'\coordinate (', name, ') {}', file=f)
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
			f.write('node ')
			ctype, images, _ = aut.orbit_type(path, intersection)
			if dashed[-1] and aut.orbit_type(path, intersection)[0].is_type_B():
				f.write('[characteristic] ')
			f.write('{' + str(leafnum) + '}')
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

# if __name__ == "__main__":
	# aut = Automorphism.from_file('../thompson/examples/semi_inf_c.aut')
	# dump_tikz_code(aut, 'output.tikz')