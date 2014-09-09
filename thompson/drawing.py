"""A collection of tools to help draw trees etc. as `SVG files <http://en.wikipedia.org/wiki/Scalable_Vector_Graphics>`_.

.. testsetup:: 
	
	from thompson.drawing import *
	import svgwrite
"""

import codecs
from collections import namedtuple
from functools import wraps
from numbers import Number
from pathlib import Path
from shutil import copyfile

import svgwrite
import pkg_resources

from .constants import *

__all__ = ['Coord', 'creates_SVG', 'new_drawing', 'set_size']
__all__.sort()

_css_cache = None
_css_file = None

class Coord(namedtuple('BaseCoord', 'x y')):
	"""A :py:func:`namedtuple <py3:collections.namedtuple>` representing a coordinate :math:`(x, y)`. Some of Python's magic methods are implemented so that points can be conveniently used in arithmetic.
	"""
	#Like vanilla namedtuples, do not allow Coord instances to have any instance attributes.
	__slots__ = ()
	
	def __new__(cls, x, y=None, scale=GRID_SPACING):
		"""Coord instances scale grid units (that trees use) up to `SVG units <http://svgwrite.readthedocs.org/en/latest/overview.html#units>`_ (that the drawing uses). By default, the scale factor is :py:data:`~thompson.constants.GRID_SPACING` (40).
		
			>>> x = Coord(1, 2)
			>>> x == (40, 80)
			True
		
		Coords can also be created from an object with ``x`` and ``y`` attributes, or from an object supporting indices 0 and 1.
		
			>>> class C: pass
			>>> c = C()
			>>> c.x, c.y = 1, 2
			>>> x == Coord(c)
			True
			>>> l = [1, 2]
			>>> x == Coord(l)
			True
		
		The scale argument is usually omitted. It is used in :meth:`__add__` and :meth:`__mul__` to make sure we don't scale by another factor of 40 when creating a new :class:`Coord`.
		"""
		if y is None:
			try:
				if hasattr(x, 'x'):
					x, y = x.x, x.y
				else:
					x, y = x[0], x[1]
			except (AttributeError, IndexError):
				raise ValueError("Second argument omitted, but first argument does not represent a coordinate (x, y).")
			
		self = super(Coord, cls).__new__(cls, x*scale, y*scale)
		return self
	
	def __add__(self, other):
		"""Coordinates can be added elementwise.
		
			>>> Coord(1, 2) + Coord(3, 4) == Coord(4, 6)
			True
		"""
		return type(self)(self.x + other.x, self.y + other.y, scale=1)
	
	def __mul__(self, other):					#self * other
		"""Coordinates can be multiplied by an integer---in fact, by any :py:class:`Number <py3:numbers.Number>`.
			
			>>> Coord(1, 2) * 10 == Coord(10, 20) == 10 * Coord(1, 2)
			True
		"""
		if isinstance(other, Number):
			return type(self)(other * self.x, other * self.y, scale=1)
		else:
			return NotImplemented
	
	def __rmul__(self, other):					#other * self
		#Ensure that other * self == self * other
		return self.__mul__(other)
	
	@classmethod
	def length(cls, x, scale=GRID_SPACING):
		"""A one-dimensional version of Coord. Scales up its argument by :data:`~thompson.constants.GRID_SPACING`."""
		#This doesn't require a new class because addition, multiplication etc work as you'd expect. But we have to tell Python how to add coordinates (tuple addition is concatenation) and scale them (tuple scaling is repeition).
		return x*scale

def creates_SVG(f):
	"""Decorator. If *f* is a function that returns a :class:`Group <svgwrite:svgwrite.container.Group>`, applying this decorator modifies *f* to optionally embed that group in an SVG file. For example:
	
		>>> @creates_SVG
		... def f(**kwargs):
		... 	return svgwrite.container.Group()
		... 
		>>> f()
		<svgwrite.container.Group object at 0x...>
		>>> dwg, canvas = f(filename='output.svg')
		>>> dwg, canvas
		(<svgwrite.drawing.Drawing object at 0x...>, <svgwrite.container.Group object at 0x...>)
		>>> dwg.filename
		'output.svg'
	
	.. testcleanup::
		
		import os
		os.remove('output.svg')
	
	When f is defined, its arguments should include ``**kwargs``.
	
	When *f* is called with a  *filename* keyword argument, :meth:`new_drawing` is called using the entirety of ``**kwargs``. (This means that other arguments can be passed :meth:`new_drawing`.) Next, *f* executes as normal, producing an SVG group. If no filename argument is present, *f* behaves exactly as if it were undecorated.
	
	This is meant to decorate functions that produce components to be used in constructing larger graphics. We encourage such functions to use :func:`set_size` before returning.
	"""
	@wraps(f)
	def output(*args, **kwargs):
		if 'filename' in kwargs:
			dwg, canvas = new_drawing(**kwargs)
		group = f(*args, **kwargs)
		if 'filename' in kwargs:
			canvas.add(group)
			dwg.save()
			return dwg, canvas
		else:
			return group
	
	if output.__doc__ is None:
		output.__doc__ = ""
	i = get_docstring_indentation(output.__doc__)
	output.__doc__ += """\n\n""" + "\t"*i + """This method is decorated by :func:`~thompson.drawing.creates_SVG`. A *filename* argument can be supplied to render the group to an SVG file."""
	return output

def get_docstring_indentation(string):
	"""
	"""
	i = 0
	if string is None:
		return i
	try:
		line = string.splitlines()[1]
	except IndexError:
		return i
	for i, char in enumerate(line):
		if char != '\t':
			break
	else: #if we didn't break from the loop (line full of tabs)
		i += 1
	return i

def new_drawing(filename=None, num_grid_lines=0, embed_css=True, debug_css=False):
	"""Creates a new SVG :class:`Drawing <svgwrite:svgwrite.drawing.Drawing>`. A :class:`Group <svgwrite:svgwrite.container.Group>` called *canvas* is created within the drawing, positioned slightly away from the origin. The grid that :class:`Coord` describes can optionally be drawn.
	
	CSS Selectors:
	
	- The canvas has id 'canvas'.
	- If drawn, the grid is contained in a group with class 'grid'.
		- Axis lines have class 'axis'.
		- Every fifth line has class 'major'.
	
	:param str filename: The filename associated with the drawing. **NB**: no file is created on disk until :py:meth:`dwg.save <svgwrite:svgwrite.drawing.Drawing.save>` is called. The *filename* may be omitted or be given as ``None``; if so, dwg.filename should be set before saving the file to disk.
	:param int num_grid_lines: Draws the given number of lines of the grid coordinate system, as well as axes.
	:param bool embed_css: If true, CSS styles are embedded in a ``<style>`` tag. If false, a ``<link>`` to an external CSS file is made. If in addition *filename* is provided, the CSS file is copied to the directory containing *filename*.
	:return: A 2-tuple (drawing, canvas).
	"""
	dwg = svgwrite.Drawing(filename)
	include_css(dwg, embed_css)
	include_markers(dwg)
	
	canvas = svgwrite.container.Group(id='canvas')
	canvas.translate(Coord(1, 1))
	dwg.add(canvas)
	
	if num_grid_lines > 0:
		grid = svgwrite.container.Group(class_='grid debug')
		for i in range(num_grid_lines+1):
			x = svgwrite.shapes.Line(Coord(0, i), Coord(num_grid_lines, i))
			y = svgwrite.shapes.Line(Coord(i, 0), Coord(i, num_grid_lines))
			if i == 0:
				x['class'] = y['class'] = 'axis'
			elif i % 5 == 0:
				x['class'] = y['class'] = 'major'
			grid.add(x)
			grid.add(y)
		canvas.add(grid)
	return dwg, canvas

def include_css(dwg, embed):
	"""Ensures that *dwg* has access to the CSS styles needed to render trees etc. Called by :func:`new_drawing`."""
	global _css_cache, _css_file
	_css_file = _css_file or pkg_resources.resource_filename('thompson', 'drawing.css')
	if embed:
		if _css_cache is None:
			with open(_css_file, 'rt', encoding='utf-8') as f:
				_css_cache = svgwrite.container.Style(f.read())
		dwg.defs.add(_css_cache)
	else:
		dwg.add_stylesheet('drawing.css')
		if dwg.filename is not None:
			dest = Path(dwg.filename)
			dest.name = 'drawing.css'
			#TODO: only copy over the css file into a directory once. Keep a set of paths it's been copied into?
			copyfile(_css_file, dest)

def include_markers(dwg):
	"""Generates some useful markers and embeds them in the given :class:`Drawing <svgwrite:svgwrite.drawing.Drawing>`."""
	#ARROWHEAD
	arrowhead = svgwrite.container.Marker(
		insert=(9, 5), size=(10, 10), orient='auto', markerUnits='strokeWidth',
		id='arrowhead'
	)
	arrowhead.add(svgwrite.path.Path("M 0 0 L 10 5 L 0 10 z"))
	dwg.defs.add(arrowhead)
	
	#TICKMARK
	tickmark = svgwrite.container.Marker(
		insert=(1, 5), size=(2, 10), orient='auto', markerUnits='strokeWidth',
		id='tickmark'
	)
	tickmark.add(svgwrite.shapes.Line((0, 0), (0, 10), stroke_width="2"))
	dwg.defs.add(tickmark)
	
	#DISCONTINUITY
	discontinuity = svgwrite.container.Marker(
		insert=(4, 4), size=(8, 8), orient='auto', markerUnits='strokeWidth',
		id='discontinuity'
	)
	discontinuity.add(svgwrite.shapes.Circle((4, 4), 3.5))
	dwg.defs.add(discontinuity)
	
	#CONTINUITY
	continuity = svgwrite.container.Marker(
		insert=(4, 4), size=(8, 8), orient='auto', markerUnits='strokeWidth',
		id='continuity'
	)
	continuity.add(svgwrite.shapes.Circle((4, 4), 3.5))
	dwg.defs.add(continuity)
	
	#CROSSHAIR
	crosshair = svgwrite.container.Marker(
		insert=(5, 5), size=(10, 10), markerUnits='strokeWidth',
		id='crosshair'
	)
	crosshair.add(svgwrite.path.Path("M 0 5 L 10 5 M 5 0 L 5 10"))
	dwg.defs.add(crosshair)

def set_size(g, size, offset=None):
	"""SVG Groups are used to store various components of drawings. For instance, a :class:`~thompson.tree_pair.TreePair` is rendered as a group containing two groups, one for each tree. Once such a group has been created, it is useful to align these to the grid of :class:`Coord` inates.
	
	This function sets the size attribute of the group to the given *size*. This attribute is not provided by :module:`svgwrite`, but is a custom addition internally used by ``thompson``. A rectangle with CSS class 'debug' is added to the group, in order to visualise its size. If provided, the group is translated by the *offset* argument, though the debug rectangle is not translated.
	
	:param size: A :class:`Coord` inate.
	:param offset: Optional: a :class:`Coord` inate.
	
	**Example.**
	
	.. figure:: examples/set_size_example.svg
		
		The effect of the *offset* argument. Left: no offset; Right: offset ``Coord(0.5, 0.5)``. [:download:`Source code <examples/set_size_example.py>`].
	"""
	if offset is not None:
		g.translate(offset)
		insert = -1*offset
	else:
		insert = (0, 0)
	r = svgwrite.shapes.Rect(insert=insert, size=size, class_='debug')
	g.elements.insert(0, r)
	g.size = size