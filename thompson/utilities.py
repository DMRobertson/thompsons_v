"""Functions which use the data structures of ``thompson'' but aren't part of the core functionality."""

	# level/.style={sibling distance=(1mm + 60mm*0.5^#1)},
head = r'''\tikzset{
	level distance=8mm,
	level 1/.style={sibling distance=30mm},
	level 2/.style={sibling distance=15mm},
	level 3/.style={sibling distance=10mm},
	level 4/.style={sibling distance=5mm},
	level 5/.style={sibling distance=3.5mm},
	characteristic/.style={red, draw, solid}
}
'''
tail = ''''''
#tail used to be '''\end{tikzpicture}''', but we removed that so tweaks can be made to the drawing from the main document, e.g. \begin{tikzpicture}[scale=0.5]

def generate_tikz_code(aut, filename, domain=None, name=''):
	r"""Saves instructions in the given *filename* for TikZ to draw the given automorphism. The leaves of the domain tree are the elements of *domain*. The arrow is labelled with the contents of *name*; this can include TeX mathematics syntax, e.g. ``name=r'\alpha'``.
	
	The difference forests domain - range and range - domain are drawn with dashed lines. Attractors and repellers are drawn in red.
	"""
	#1. Compute X = L(domain \intersect range)
	if domain is None:
		domain = aut.domain
	range = aut.image_of_set(domain)
	X = domain.copy()
	X += range
	finished = False
	while not finished:
		i, j = X.test_free()
		if i == j == -1:
			finished = True
		else:
			del X[j] #this is never above X[i], see the docstring for test_free
	X.sort()
	
	#2. Generate the tikz code
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
	depth = 1
	dashed = []
	
	print(r'\begin{scope}[', extra, ']', file=f)
	print(r'\coordinate (', name, ') {}', file=f)
	
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