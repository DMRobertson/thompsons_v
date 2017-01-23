from ..generators import Generators
from ..utilities  import intersection_of_trees, intersection_from_domain, handle_domain
import pkg_resources

#TODO load these when required, not at import time
graph_header      = pkg_resources.resource_string( __name__, 'graph_head.tex').decode('utf-8')
graph_footer      = pkg_resources.resource_string( __name__, 'graph_foot.tex').decode('utf-8')
standalone_header = pkg_resources.resource_string( __name__, 'standalone_head.tex').decode('utf-8')
standalone_footer = pkg_resources.resource_string( __name__, 'standalone_foot.tex').decode('utf-8')
diagram_header    = pkg_resources.resource_string( __name__, 'diagram_head.tex').decode('utf-8')
diagram_footer    = pkg_resources.resource_string( __name__, 'diagram_foot.tex').decode('utf-8')

def write_tikz_code(aut,
	domain='wrt QNB',
	dest='tikz_generated.tex',
	horiz=True,
	name='',
	standalone=True,
	draw_revealing=True
):
	#1. Decide which domain to use for plotting.
	domain = handle_domain(domain, aut)
	range, intersection = intersection_from_domain(domain, aut)
	# Number the leaves.
	domain = [ (w, i + 1) for (i, w) in enumerate(domain) ]
	range  = [ (w, i + 1) for (i, w) in enumerate(range)  ]
	#Order the range using Higman's words
	range.sort()

	with open(dest, 'wt', encoding='utf-8') as f:
		if standalone:
			f.write(standalone_header)
		f.write(diagram_header)

		write_basis(f, aut, domain, True, horiz, intersection, draw_revealing)
		write_basis(f, aut, range, False, horiz, intersection, draw_revealing)

		write_arrow(f, horiz, name)
		f.write(diagram_footer)
		if standalone:
			f.write(standalone_footer)

def write_arrow(f, horiz, name):
	"""The connecting arrow from the domain forest to range"""
	if name.startswith('$') and not name.endswith('$'):
		raise ValueError("Arrow names must end with a $ if they begin with a $.")

	f.write("\\draw[->, thick, shorten >=0.5em, shorten <=0.5em]\n")
	if horiz:
		f.write("\tlet \\p1=(domain.east), \\p2=(range.west), \\n1={max(\\y1,\\y2)} in\n")
		f.write("\t\t(\\x1, \\n1) -- node[auto] {{{}}} (\\x2, \\n1);\n".format(name))
	else:
		f.write("\t(domain) -- node[auto] {{{}}} (range);\n".format(name))

def name(word, for_domain, options):
	"""Introduce a node with the given key/value options."""
	prefix = 'd' if for_domain else 'r'
	label = prefix + str(word).replace(" ", "")
	pairs = []
	for key, value in options.items():
		if value is None:
			pairs.append(key)
		else:
			pairs.append('{}={}'.format(key, value))
	if pairs:
		label += " [{}]".format( ', '.join(pairs) )
	return label

def write_basis(f, aut, basis, for_domain, horiz, intersection, draw_revealing):
	write_enclosing_node(f, for_domain, horiz)
	f.write(graph_header)

	for i, (word, label) in enumerate(basis):
		if draw_revealing:
			highlight = is_repatt(word, intersection, for_domain, aut)
		else:
			highlight = False
		write_word(f, word, label, for_domain, intersection, highlight)

	f.write(graph_footer)
	f.write("};\n")

def write_word(f, word, label, for_domain, intersection, highlight):
	below_intersection = False
	for subword in word.subwords():
		options = {}
		f.write("\t\t\t")

		if len(subword) == 1:
			options["root"] = None
		else:
			f.write(" -- ")

		if below_intersection:
			f.write("[component] ")
			if highlight:
				#TODO I'm using two different ways to style an edge here---doesn't make sense.
				options["target edge style"] = "spine"

		if subword == word:
			options["leaf"] = None
			options["label"] = "below:" + str(label)
			if highlight:
				options["label"] = "{[repatt label]" + options["label"] + "}"
				options["repatt"] = None

		f.write(name(subword, for_domain, options))

		if not below_intersection:
			below_intersection = subword in intersection

		if subword == word:
			f.write(",")
		f.write("\n")

def write_enclosing_node(f, for_domain, horiz):
	"""We draw each forest as a tikzpicture inside a node."""
	f.write("\\node (" + ('domain' if for_domain else 'range') +  ") ")
	if not for_domain:
		if horiz:
			f.write('[right=5em of domain.north east, anchor=north west] ')
		else:
			f.write('[below=5em of domain.south, anchor=north]')
	f.write("{\n")

def is_repatt(word, intersection, for_domain, aut):
	if word in intersection:
		return False
	ctype = aut.orbit_type(word, intersection)[0]
	if not ctype.is_type_B():
		return False
	if (ctype.characteristic[0] < 0) != for_domain:
		return False
	#TODO check that the root you stop at is above where you started
	return True
