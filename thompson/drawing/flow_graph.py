from ..orbits import ComponentType

combining_overline = "\u0305"

def finite_orbits(aut):
	for _, orbits in aut.periodic_orbits.items():
		yield from orbits

def name(word):
	return str(word).replace(" ", "")

def flow_graph_code(aut, basis=None, filename="flow_graph.dot"):
	"""Generate DOT code for the flow graph for *aut* with respect to the quasinormal basis."""
	with open(filename, "wt") as f:
		f.write("digraph {")
		
		#Periodic bits, one line per orbit
		for orbit in finite_orbits(aut):
			print("orbit: ", orbit)
			f.write("\n\t")
			for word in orbit:
				f.write( name(word) + " -> " )
			f.write( name(orbit[0]) + "\n" )
		
		#Reg infinite bits
		init, term = aut.semi_infinite_end_points()
		f.write("\t#Initial\n")
		for word in init:
			if word not in aut.quasinormal_basis:
				continue
			f.write("\t" + name(word) + "\n" )
		f.write("\t#Terminals\n")
		for word in term:
			if word not in aut.quasinormal_basis:
				continue
			f.write("\t" + name(word) + "\n" )
		f.write("\t#Flow lines\n")
		
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
			
			root, tail = term.test_above(d)
			f.write("\t" + name(root))
			index = 1
			while index + 2 in images:
				f.write( " -> " + name(images[index]) )
				index += 1
			element = images[index]
			root, tail = init.test_above(element)
			f.write(" -> " + name(root))
			f.write("\n")
			
		
		
		
		
		#End of graph
		f.write("}\n")
