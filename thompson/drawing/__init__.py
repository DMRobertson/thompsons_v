from .plot     import plot as plot_svg
from .treepair import write_tikz_code

from subprocess import call, check_call
from tempfile   import mkstemp, mkdtemp

import os
import sys

def display_file(filepath):
	"""From http://stackoverflow.com/a/435669. Opens the given file with the OS's default application."""
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
		display_file(dest)
	return dest

def forest(aut, jobname=None, name='', display=True, horiz=True):
	outdir = mkdtemp()
	if jobname is None:
		jobname = 'forest_diagram'
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
		display_file( pdf )
	return pdf
	
