"""Functions which use the data structures of ``thompson'' but aren't part of the core functionality."""

head_sc = r'''\documentclass{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}
'''

tail_sc = r'''\end{tikzpicture}
\end{document}
'''

	# level/.style={sibling distance=(1mm + 60mm*0.5^#1)},
head = r'''\tikzset{
	level distance=8mm,
	level 1/.style={sibling distance=35mm},
	level 2/.style={sibling distance=20mm},
	level 3/.style={sibling distance=12mm},
	level 4/.style={sibling distance=6mm},
	level 5/.style={sibling distance=3.5mm},
	characteristic/.style={red, draw, solid}
}
'''
tail = ''''''
#tail used to be '''\end{tikzpicture}''', but we removed that so tweaks can be made to the drawing from the main document, e.g. \begin{tikzpicture}[scale=0.5]

def generate_tikz_code(aut, filename, domain=None, name='', self_contained=False):
	r""".. caution:: This is an experimental feature based partially on [SD10]_.
	
	Saves instructions in the given *filename* for TikZ to draw the given automorphism. The leaves of the domain tree are the elements of *domain*. The arrow is labelled with the contents of *name*; this can include TeX mathematics syntax, e.g. ``name=r'$\alpha$'``.
	
	The difference forests domain - range and range - domain are drawn with dashed lines. Attractors and repellers are drawn in red.
	"""
	#1. Compute X = L(domain \intersect range)
	domain = handle_domain(domain, aut)
	range, X = intersection_from_domain(domain, aut)
	range = aut.image_of_set(domain)
	
	#2. Generate the tikz code
	with open(filename, 'wt') as f:
		if self_contained:
			f.write(head_sc)
		f.write(head)
		basis_to_tree(aut, domain, range, X, f, name='domain')
		
		y = max(len(w) for w in domain) - 1
		z = max(len(w) for w in domain) - 1
		depth = -min(y, z)/5
		coordend = r'\textwidth, ' + str(depth*8) + 'mm)'
		print(r'\draw[->] (0.25' + coordend, '-- node[auto]{', name, '} (0.35' + coordend, ';\n', file=f)
		
		basis_to_tree(aut, range, domain, X, f, name='range', extra=r'xshift=0.6\textwidth')
		f.write(tail)
		if self_contained:
			f.write(tail_sc)

from subprocess   import call, check_call
from tempfile     import mktemp, mkdtemp

def write_tikz_code(self, filename, domain='wrt QNB', name='', self_contained=False):
	generate_tikz_code(self, filename, domain, name, self_contained)
write_tikz_code.__doc__ = generate_tikz_code.__doc__
	
def render(self, jobname=None, domain='wrt QNB', name='', display=True):
	"""Uses :meth:`write_tikz_code` and a call to ``pdflatex`` to generate a PDF drawing of the given automorphism. If the *display* argument is ``True``, then the system will attempt to open the PDF using the OS's default application for PDFs.
	
	.. caution:: This is an experimental feature based on [SD10]_.
	"""
	outdir = mkdtemp()
	if jobname is None:
		jobname = 'tikz_code'
	tex_file = os.path.join(outdir, jobname + '.tex')
	pdf_file = os.path.join(outdir, jobname + '.pdf')
	
	self.write_tikz_code(tex_file, domain=domain, name=name, self_contained=True)
	check_call(['pdflatex',
		'-output-directory=' + outdir,
		'-interaction=batchmode',
		'-no-shell-escape',
		 tex_file
	])
	if display:
		display_file(pdf_file)
	return pdf_file
