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
	return range, aut._intersection_of_trees(domain, range)



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

def is_repeller_attractor(path, aut, in_difference, other, intersection, in_domain):
	if not in_difference:
		return False
	ctype = aut.orbit_type(path, intersection)[0]
	if not ctype.is_type_B():
		return False
	if not (ctype.characteristic.power < 0) == in_domain:
		return False
	ancestor = aut.repeated_image(path, -ctype.characteristic.power)
	assert ancestor.is_above(path)
	return ancestor in other
