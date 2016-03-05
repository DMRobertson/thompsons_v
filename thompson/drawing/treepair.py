from ..generators import Generators
from ..utilities  import intersection_of_trees
import pkg_resources

graph_header      = pkg_resources.resource_string( __name__, 'graph_head.tex').decode('utf-8')
graph_footer      = pkg_resources.resource_string( __name__, 'graph_foot.tex').decode('utf-8')
standalone_header = pkg_resources.resource_string( __name__, 'standalone_head.tex').decode('utf-8')
standalone_footer = pkg_resources.resource_string( __name__, 'standalone_foot.tex').decode('utf-8')
diagram_header = pkg_resources.resource_string( __name__, 'diagram_head.tex').decode('utf-8')
diagram_footer = pkg_resources.resource_string( __name__, 'diagram_foot.tex').decode('utf-8')

#todo highlight repellers and attractors

def write_tikz_code(aut,
	dest='tikz_generated.tex',
	horiz=True,
	name='',
	standalone=True, 		
	annotate=True
):
	with open(dest, 'wt', encoding='utf-8') as f:
		if standalone:
			f.write(standalone_header)
		f.write(diagram_header)
		
		if annotate:
			intersect = intersection_of_trees(aut.domain, aut.range)
			compare = (intersect, aut)
		else:
			compare = None
		
		domain = [ (w, i + 1) for (i, w) in enumerate(aut.domain) ]
		write_basis(f, domain, True, horiz, compare)
		
		range  = [ (w, i + 1) for (i, w) in enumerate(aut.range)  ]
		range.sort()
		write_basis(f, range, False, horiz, compare)
		
		write_arrow(f, horiz, name)
		f.write(diagram_footer)
		if standalone:
			f.write(standalone_footer)

def write_arrow(f, horiz, name):
	f.write("\draw[->, thick, shorten >=0.5em, shorten <=0.5em]\n")
	if horiz:
		f.write("\tlet \\p1=(domain.east), \\p2=(range.west), \\n1={max(\y1,\y2)} in\n")
		f.write("\t\t(\\x1, \\n1) -- node[auto] {{{}}} (\\x2, \\n1);\n".format(name))
	else:
		f.write("\t(domain) -- node[auto] {{{}}} (range);\n".format(name))
		
def name(word, domain, options):
	prefix = 'd' if domain else 'r'
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

def pairs(iterator):
	this = next(iterator)
	while True:
		try:
			that = next(iterator)
		except StopIteration:
			yield (this, None)
			break
		else:
			yield (this, that)
			this = that

def write_basis(f, basis, domain, horiz, compare=None):
	if compare is None:
		draw_revealing = False
	else:
		draw_revealing = True
		intersection, aut = compare
			
	write_enclosing_node(f, domain, horiz)
	f.write(graph_header)
	
	repatts = []
	
	for i, (word, label) in enumerate(basis):
		highlight = is_repatt(word, intersection, domain, aut)
		write_word(f, word, label, draw_revealing, highlight, domain, intersection)

	f.write(graph_footer)
	f.write("};\n")

def write_word(f, word, label, draw_revealing, highlight, domain, intersection):
	below_intersection = False
	for subword in word.subwords():
		options = {}
		f.write("\t\t\t")
		
		if len(subword) == 1:
			options["root"] = None
		else:
			f.write(" -- ")
		
		if draw_revealing and below_intersection:
			f.write("[component] ")
			if highlight:
				options["target edge style"] = "spine"
		
		if subword == word:
			options["leaf"] = None
			options["label"] = "below:" + str(label)
			if draw_revealing and highlight:
				options["label"] = "{[repatt label]" + options["label"] + "}"
				options["repatt"] = None
		
		f.write(name(subword, domain, options))
		
		if draw_revealing and not below_intersection:
			below_intersection = subword in intersection
		
		if subword == word:
			f.write(",")
		f.write("\n")

def write_enclosing_node(f, domain, horiz):
	f.write("\\node (" + ('domain' if domain else 'range') +  ") ")
	if not domain:
		if horiz:
			f.write('[right=5em of domain.north east, anchor=north west] ')
		else:
			f.write('[below=5em of domain.south, anchor=north]')		
	f.write("{\n")

def is_repatt(word, intersection, domain, aut):
	if word in intersection:
		return False
	ctype = aut.orbit_type(word, intersection)[0]
	return ctype.is_type_B() and (ctype.characteristic[0] < 0) == domain
