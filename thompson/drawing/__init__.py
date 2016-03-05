from .plot     import plot as plot_svg
from .treepair import write_tikz_code

from subprocess import call, check_call
from tempfile   import mkstemp, mkdtemp

import os
import sys

def display_file(filepath, format=None):
	if in_ipynb():
		from IPython.display import Image, SVG
		if format == 'svg':
			return SVG(filename=filepath)
		elif format == 'pdf':
			png_file = filepath[:-4] + '.png'
			check_call(['convert',
				'-density', '180',
				filepath,
				'-quality','90',
				png_file
			])
			return Image(filename=png_file)
		else:
			raise NotImplementedError
	
	"""From http://stackoverflow.com/a/435669.
	Opens the given file with the OS's default application."""
	if sys.platform.startswith('darwin'):
		call(('open', filepath))
	elif os.name == 'nt':
		os.startfile(filepath)
	elif os.name == 'posix':
		call(('xdg-open', filepath))
	else:
		raise NotImplementedError

def plot(aut, dest=None, display=True, endpoints=False):
	if dest is None:
		#Write to a temporary file
		dest = mkstemp()[1]
	plot_svg(aut, dest, endpoints)
	
	if display:
		return display_file(dest, format='svg')
	return dest

def forest(aut, jobname=None, name='', display=True, horiz=True):
	if jobname is None:
		outdir = mkdtemp()
		jobname = 'forest_diagram'
	else:
		outdir = '.'
	tex = os.path.join(outdir, jobname + '.tex')
	pdf = os.path.join(outdir, jobname + '.pdf')
	
	write_tikz_code(aut, tex, horiz, name)
	check_call(['lualatex',
		'-output-directory=' + outdir,
		'-interaction=batchmode',
		'-no-shell-escape',
		 tex
	])
	
	if display:
		return display_file(pdf, format='pdf')
	return pdf

def in_ipynb():
	"""From a comment http://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook#comment53916407_24937408"""
	try:
		ipy = get_ipython()
	except NameError:
		return False
	else:
		return type(ipy).__name__ == 'ZMQInteractiveShell'
		
