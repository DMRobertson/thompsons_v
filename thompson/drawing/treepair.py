from thompson.examples import random_automorphism
import pkg_resources

def write_tikz_code(aut, dest='tikz_generated.tex'):
	header = pkg_resources.resource_filename('thompson.drawing', 'tikz_template_head.tex')
	footer = pkg_resources.resource_filename('thompson.drawing', 'tikz_template_foot.tex')
	with open(dest, 'wt') as f:
		append_to(header, f)
		
		for i in range(0, aut.signature.alphabet_size):
			f.write( "\tdx{} [root]\n".format( i+1 ) )
			f.write( "\trx{} [root]\n".format( i+1 ) )
		
		write_basis(f, aut.domain, True)
		write_basis(f, aut.range, False)
		
		

				
		append_to(footer, f)

def append_to(src, dest):
	with open(src, 'rt') as f:
		for line in f:
			dest.write(line)

def name(word, domain):
	prefix = 'd' if domain else 'r'
	return prefix + str(word).replace(" ", "")

def write_basis(f, basis, domain):
	for i, word in enumerate(basis):
		for subword in word.subwords():
			node_name = name(subword, domain=True)
			first = len(subword) == 1
			last  = len(subword) == len(word)
			f.write("\t")
			f.write(node_name)
			if len(subword) > 1:
				index = -subword[-1]
				f.write( " [desired child index={}] ".format(index) )
			if last:
				f.write( " [leaf, as=${}$],\n".format(i+1) )

# this almost works
"""		
		for i, word in enumerate(sorted(aut.range)):
			names = [name(subword, domain=False) for subword in word.subwords()]
			names[-1] += " [leaf, as=${}$],\n".format(i+1)
			f.write("\t" + " -- ".join(names))
"""

if __name__ == "__main__":
	f = random_automorphism()
	write_tikz_code(f)
