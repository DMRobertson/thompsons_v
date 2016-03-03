from thompson import Generators
from thompson.examples import random_automorphism
import pkg_resources

graph_header      = pkg_resources.resource_string( __name__, 'graph_head.tex').decode('utf-8')
graph_footer      = pkg_resources.resource_string( __name__, 'graph_foot.tex').decode('utf-8')
standalone_header = pkg_resources.resource_string( __name__, 'standalone_head.tex').decode('utf-8')
standalone_footer = pkg_resources.resource_string( __name__, 'standalone_foot.tex').decode('utf-8')
diagram_header = pkg_resources.resource_string( __name__, 'diagram_head.tex').decode('utf-8')
diagram_footer = pkg_resources.resource_string( __name__, 'diagram_foot.tex').decode('utf-8')

#todo highlight repellers and attractors

def write_tikz_code(aut, dest='tikz_generated.tex', horiz=True, name='', standalone=True):
	with open(dest, 'wt', encoding='utf-8') as f:
		if standalone:
			f.write(standalone_header)
		f.write(diagram_header)
		
		domain = [ (w, i + 1) for (i, w) in enumerate(aut.domain) ]
		write_basis(f, domain, True, horiz)
		
		range  = [ (w, i + 1) for (i, w) in enumerate(aut.range)  ]
		range.sort()
		write_basis(f, range, False, horiz)
		
		f.write("\draw[->, thick, shorten >=0.5em, shorten <=0.5em]\n")
		if horiz:
			f.write("\tlet \\p1=(domain.east), \\p2=(range.west), \\n1={max(\y1,\y2)} in\n")
			f.write("\t\t(\\x1, \\n1) -- node[auto] {{{}}} (\\x2, \\n1);\n".format(name))
		else:
			f.write("\t(domain) -- node[auto] {{{}}} (range);\n".format(name))
		
		f.write(diagram_footer)
		if standalone:
			f.write(standalone_footer)

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

def write_basis(f, basis, domain, horiz):
	f.write("\\node (" + ('domain' if domain else 'range') +  ") ")
	if not domain:
		if horiz:
			f.write('[right=5em of domain.north east, anchor=north west] ')
		else:
			f.write('[below=5em of domain.south, anchor=north]')
	
	f.write("{\n")
	f.write(graph_header)
	
	placed = set()
	for i, (word, label) in enumerate(basis):
		first = True
		for parent, child in pairs(word.subwords()):
			if child in placed:
				continue
			options = {}
			
			f.write("\t\t\t")
			if not first:
				f.write("\t-- ")
			else:
				first = False
			
			if len(parent) == 1:
				options["root"] = None
			if parent == word:
				options["leaf"] = None
				options["label"] = "below:" + str(label)
			
			f.write(name(parent, domain, options))
			placed.add(parent)
			
			if parent == word:
				if i == len(basis) - 1:
					f.write(";")
				else:
					f.write(",")
			f.write("\n")
	f.write(graph_footer)
	f.write("};\n")

if __name__ == "__main__":
	f = random_automorphism()
	print(f)
	write_tikz_code(f, horiz=True, name='$f$')
