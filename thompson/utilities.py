from .generators import Generators

"""Functions whose definitions we dont want cluttering other modules' source code."""

def intersection_of_trees(domain, range):
	r"""Given the leaves A and B of two trees (i.e. bases), computes the leaves (basis) X of the tree intersection A & B.
	"""
	assert domain.signature == range.signature
	basis = Generators.standard_basis(domain.signature)
	i = 0
	while i < len(basis):
		b = basis[i]
		if b in domain or b in range:
			i += 1
		else:
			basis.expand(i)
	return basis

def handle_domain(domain, aut):
	if isinstance(domain, str):
		d = domain.lower().strip()
		if d == 'wrt qnb':
			domain = aut.quasinormal_basis.minimal_expansion_for(aut)
		elif d == 'minimal':
			domain = aut.domain
		else:
			raise ValueError("Don't know what {!r} means when specifying a domain".format(domain))
	return domain

def intersection_from_domain(domain, aut):
	"""Given a basis :math:`X` and an automorphism :math:`f`, we can always construct the minimal expansion :math:`Y` of :math:`X`. This function provides an inverse of sorts: given :math:`Y`, let :math:`Z = f(Y)` and define :math:`X` to be the basis corresponding to the intersection :math:`Y \cap Z` of trees."""
	range = aut.image_of_set(domain)
	return range, intersection_of_trees(domain, range)


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
			if is_repeller_attractor(path, aut, dashed[-1], other, intersection, extra == ''):
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


