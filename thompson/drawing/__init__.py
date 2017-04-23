"""
.. testsetup::

	from thompson.drawing import *
"""

from .plot       import plot as plot_svg
from .forest     import forest_code
from .flow_graph import flow_graph_code

from subprocess import call, check_call
from tempfile   import mkstemp, mkdtemp

from io import StringIO
import os
import shutil
import sys

__all__ = ["display_file", "forest", "forest_code", "plot", "flow"]

def display_file(filepath, format=None, scale=1.0, verbose=False):
	"""Display the image at *filepath* to the user. This function behaves differently, depending on whether or not we execute it in a Jupyter notebook.

	If we are **not** in a notebook, this function opens the given image using the operating system's default application for that file.

	If we **are** in a notebook, this returns an IPython object corresponding to the given image.
	If this object is the last expression of a notebook code block, the image will be displayed.
	The image is handled differently depending on its *format*, which must be specified when the function is called in a notebook.
	Only the formats `'svg'` and `'pdf'` are accepted.

	.. note::

		PDF files are displayed in a notebook as a rendered PNG. The conversion is made using the ``convert`` program provided by `ImageMagick <http://www.imagemagick.org/script/index.php>`_, which must be available on the `PATH <https://en.wikipedia.org/wiki/PATH_(variable)>`_ for this function to work with PDF files. The optional *scale* argument can be used to control how big the rendered PNG is.

	.. todo:: use a Python binding to ImageMagick rather than just shelling out?
	"""
	if in_ipynb():
		print(filepath)
		from IPython.display import Image, SVG
		format = format.lower()
		if format == 'svg':
			return SVG(filename=filepath)
		elif format == 'png':
			return Image(filename=filepath)
		elif format == 'pdf':
			png_file = filepath[:-4] + '.png'
			check_call(['convert',
				'-density', str(100*scale),
				filepath,
				'-quality','90',
				png_file
			])
			return Image(filename=png_file)
		else:
			raise NotImplementedError( "Don't know how to display format {}".format(format) )

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
	return filepath

def plot(*auts, dest=None, display=True, diagonal=False, endpoints=False):
	r"""Plots the given :class:`automorphism <thompson.automorphism.Automorphism>` s as a function :math:`[0, 1] \to [0, 1]`. The image is rendered as an SVG using `svgwrite`.

	:param str dest: the destination filepath to save the SVG to. If `None`, the SVG is saved to a temporary file location.
	:param bool display: if True, automatically call :func:`display_file` to display the SVG to the user. Otherwise does nothing.
	:param bool diagonal: if True, draws the diagonal to highlight fixed points of *aut*.
	:param bool endpoints: if True, open and closed endpoints are drawn on the plot to emphasise any discontinuities that *aut* has.

	:returns: the filepath where the SVG was saved. If *dest* is `None` this is a temporary file; otherwise the return value is simply *dest*.
	"""
	if dest is None:
		#Write to a temporary file
		dest = mkstemp()[1]
	plot_svg(*auts, filename=dest, diagonal=diagonal, endpoints=endpoints)

	if display:
		return display_file(dest, format='svg')
	return dest

def forest(aut, jobname=None, display=True, scale=1, **kwargs):
	r"""Draws the given :class:`~thompson.automorphism.Automorphism` as a forest-pair diagram.
	The image is rendered as a PDF using the `tikz` graph drawing libraries and `lualatex`.

	:param str jobname: the destination filepath to save the PDF to. A file extension should **not** be provided. If `None`, the PDF is saved to a temporary file location.
	:param bool display: if True, automatically call :func:`display_file` to display the PDF to the user. Otherwise does nothing.
	:param float scale: In a Jupyter notebook, this controls the size of the rendered PNG image. See the note in :func:`display_file`.

	:returns: the filepath where the PDF was saved. If *dest* is `None` this is a temporary file; otherwise the return value is simply `jobname + '.pdf'`.

	.. note::

		The graph drawing is done via a TikZ and LaTeX. The source file is compiled using `lualatex`, which must be available on the `PATH <https://en.wikipedia.org/wiki/PATH_(variable)>`_ for this function to work.
	"""
	outdir = mkdtemp()
	
	specific_location = jobname is not None
	if not specific_location:
		jobname = 'forest_diagram'
	else:
		destdir = os.path.dirname(jobname)
		jobname = os.path.basename(jobname)
	tex = os.path.join(outdir, jobname + '.tex')
	pdf = os.path.join(outdir, jobname + '.pdf')
	
	with open(tex, 'wt') as texfile:
		texfile.write(forest_code(aut, **kwargs))
	check_call(['lualatex',
		'-output-directory=' + outdir,
		'-interaction=batchmode',
		'-no-shell-escape',
		 tex
	])
	
	if specific_location:
		new_pdf = os.path.join(destdir, jobname + '.pdf')
		shutil.copy(pdf, new_pdf)
		pdf = new_pdf
	
	if display:
		return display_file(pdf, format='pdf', scale=scale)
	return pdf

def forest_plain(aut, **overrides):
	"""The same as :func:`forest_code`, but this is a custom version for my own use. It passes in arguments to produce minimal TikZ code only, and draws RTL."""
	options = dict(
		include_styles = False,
		standalone = False,
		draw_revealing = False,
		LTR = False
	)
	options.update(overrides)
	return forest_code(aut, **options)

def in_ipynb():
	"""From a comment on Stack overflow:
	http://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook#comment53916407_24937408"""
	try:
		ipy = get_ipython()
	except NameError:
		return False
	else:
		return type(ipy).__name__ == 'ZMQInteractiveShell'

def flow(aut, jobname=None, display=True):
	if jobname is None:
		outdir = mkdtemp()
		jobname = 'flow_graph'
	else:
		outdir = os.path.basename(jobname)
	dot = os.path.join(outdir, jobname + ".dot")
	png = os.path.join(outdir, jobname + ".png")
	flow_graph_code(aut, dot)
	check_call([
		'dot',
		'-T', 'png',
		'-o', png,
		dot
	])
	if display:
		return display_file(png, format='png')
	return png
	
