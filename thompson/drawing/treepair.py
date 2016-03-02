from thompson import Generators
from thompson.examples import random_automorphism
import pkg_resources

def write_tikz_code(aut, dest='tikz_generated.tex', horiz=True):
	header = pkg_resources.resource_filename('thompson.drawing', 'tikz_template_head.tex')
	footer = pkg_resources.resource_filename('thompson.drawing', 'tikz_template_foot.tex')
	with open(dest, 'wt') as f:
		append_to(header, f)
		f.write( "[component direction={}]\n".format( "right" if horiz else "down" ) )
		f.write("\graph [nodes=caret] {\n")
		
		domain = [ (w, i + 1) for (i, w) in enumerate(aut.domain) ]
		write_basis(f, domain, True)
		f.write("\n")
		
		range  = [ (w, i + 1) for (i, w) in enumerate(aut.range)  ]
		range.sort()
		write_basis(f, range, False)
		
		append_to(footer, f)

def append_to(src, dest):
	with open(src, 'rt') as f:
		for line in f:
			dest.write(line)

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

def write_basis(f, basis, domain):
	subgraph_start = ( "{0} [name={0}]".format("domain" if domain else "range") )
	subgraph_start += " // [tree layout, grow=down, component direction=right, component sep=2.5em] {\n"
	f.write(subgraph_start)
	
	placed = set()
	for i, (word, label) in enumerate(basis):
		first = True
		for parent, child in pairs(word.subwords()):
			if child in placed:
				continue
			options = {}
			
			f.write("\t")
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
	
	f.write("};\n")

if __name__ == "__main__":
	f = random_automorphism()
	print(f)
	write_tikz_code(f)
