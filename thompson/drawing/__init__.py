"""
.. testsetup::
	
	from thompson.drawing import *
"""

from .plot   import plot as plot_svg
from .forest import write_tikz_code

from subprocess import call, check_call
from tempfile   import mkstemp, mkdtemp

import os
import sys

__all__ = ["display_file", "forest", "plot"]

def display_file(filepath, format=None):
	"""Display the image at *filepath* to the user. This function behaves differently, depending on whether or not we execute it in a Jupyter notebook.
	
	If we are **not** in a notebook, this function opens the given image using the operating system's default application for that file.
	
	If we **are** in a notebook, this returns an IPython object corresponding to the given image.
	If this object is the last expression of a notebook code block, the image will be displayed.
	The image is handled differently depending on its *format*, which must be specified when the function is called in a notebook.
	Only the formats `'svg'` and `'pdf'` are accepted.
	
	.. note::
		
		PDF files are displayed in a notebook as a rendered PNG. The conversion is made using the ``convert`` program provided by `ImageMagick <http://www.imagemagick.org/script/index.php>`_, which must be available on the `PATH <https://en.wikipedia.org/wiki/PATH_(variable)>`_ for this function to work with PDF files.
	
	.. todo:: use a Python binding to ImageMagick rather than just shelling out?
	"""
	if in_ipynb():
		from IPython.display import Image, SVG
		format = format.lower()
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
	r"""Plots the given :class:`automorphism <thompson.automorphism.Automorphism>` as a function :math:`[0, 1] \to [0, 1]`. The image is rendered as an SVG using `svgwrite`.
	
	:param str dest: the destination filepath to save the SVG to. If `None`, the SVG is saved to a temporary file location.
	:param bool display: if True, automatically call :func:`display_file` to display the SVG to the user. Otherwise does nothing.
	:param bool endpoints: if True, open and closed endpoints are drawn on the plot to emphasise any discontinuities that *aut* has.
	
	:returns: the filepath where the SVG was saved. If *dest* is `None` this is a temporary file; otherwise the return value is simply *dest*.
	"""
	if dest is None:
		#Write to a temporary file
		dest = mkstemp()[1]
	plot_svg(aut, dest, endpoints)
	
	if display:
		return display_file(dest, format='svg')
	return dest

def forest(aut, jobname=None, name='', display=True, horiz=True, domain=None):
	r"""Draws the given :class:`~thompson.automorphism.Automorphism` as a forest-pair diagram.
	The image is rendered as a PDF using the `tikz` graph drawing libraries and `lualatex`.
	
	:param str jobname: the destination filepath to save the PDF to. A file extension should **not** be provided. If `None`, the SVG is saved to a temporary file location.
	:param str name: The label used for the arrow between domain and range forests.
	:param bool display: if True, automatically call :func:`display_file` to display the PDF to the user. Otherwise does nothing.
	:param bool horiz: if True, draws the range forest to the right of the domain forest. If false, draws the range forest below the range forest.
	:param `~thompson.generators.Generators` domain: By default, we use the :meth:`minimal expansion <thompson.generators.Generators.minimal_expansion_for>` of the :meth:`quasi-normal basis <thompson.automorphism.Automorphism.compute_quasinormal_basis>` as the leaves of the domain forest. This can be overridden by providing a *domain* argument.
	
	:returns: the filepath where the PDF was saved. If *dest* is `None` this is a temporary file; otherwise the return value is simply `jobname + '.pdf'`.
	
	.. note::
		
		The graph drawing is done via a TikZ and LaTeX. The source file is compiled using `lualatex`, which must be available on the `PATH <https://en.wikipedia.org/wiki/PATH_(variable)>`_ for this function to work.
	"""
	if domain is None:
		domain = 'wrt QNB'
	if jobname is None:
		outdir = mkdtemp()
		jobname = 'forest_diagram'
	else:
		outdir = os.path.basename(jobname)
	tex = os.path.join(outdir, jobname + '.tex')
	pdf = os.path.join(outdir, jobname + '.pdf')
	
	write_tikz_code(aut, domain, tex, horiz, name)
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
	"""From a comment on Stack overflow:
	http://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook#comment53916407_24937408"""
	try:
		ipy = get_ipython()
	except NameError:
		return False
	else:
		return type(ipy).__name__ == 'ZMQInteractiveShell'
		
