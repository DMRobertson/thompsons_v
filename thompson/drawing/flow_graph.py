from ..orbits import ComponentType

combining_overline = "\u0305"

def finite_orbits(aut):
	for _, orbits in aut.periodic_orbits.items():
		yield from orbits

def name(word):
	return str(word).replace(" ", "")

def get_chains(aut):
	chains = []
	init, term = aut.semi_infinite_end_points()	
	#Get the sources
	domain = aut.quasinormal_basis.minimal_expansion_for(aut)
	for d in domain:
		ctype, images, _ =  aut.orbit_type(d)
		if ctype != ComponentType.complete_infinite():
			continue
		try:
			if images[-1] in domain:
				continue
		except KeyError:
			pass
		chain = []
		root, tail = term.test_above(d)
		chain.append(root)
		index = 1
		while index + 2 in images:
			chain.append(images[index])
			index += 1
		element = images[index]
		root, tail = init.test_above(element)
		chain.append(root)
		chains.append(chain)
	return chains

def flow_graph_code(aut, filename="flow_graph.dot", periodic=True):
	template = setup()
	init, term = aut.semi_infinite_end_points()
	cycles = finite_orbits(aut) if periodic else []
	with open(filename, "wt") as f:
		document = template.render(
			finite_orbits   = cycles,
			name            = name,
			attractor_roots = (i for i in init if i in aut.quasinormal_basis),
			repeller_roots  = (t for t in term if t in aut.quasinormal_basis),
			chains = get_chains(aut)
		);
		f.write(document)
	
#First time setup
template = None
def setup():
	from jinja2 import Environment, PackageLoader
	global template
	if template is not None:
		return template
	env = Environment(
		loader=PackageLoader('thompson', 'drawing'),
		lstrip_blocks = True,
		trim_blocks   = True,
	)
	template = env.get_template('flow_graph.tpl')
	return template

