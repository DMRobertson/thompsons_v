"""A collection of tools to help draw trees etc. as `SVG files <http://en.wikipedia.org/wiki/Scalable_Vector_Graphics>`_.

.. testsetup:: 
	
	from thompson.drawing import *
"""

import svgwrite
from .constants import *
from numbers import Number
from collections import namedtuple

__all__ = ['Coord', 'new_drawing', 'offset_group']
__all__.sort()

class Coord(namedtuple('BaseCoord', 'x y')):
	"""A :py:func:`namedtuple <py3:collections.namedtuple>` representing a coordinate :math:`(x, y)`. Some of Python's magic methods are implemented so that points can be conveniently used in arithmetic.
	"""
	#Like vanilla namedtuples, do not allow Coord instances to have any instance attributes.
	__slots__ = ()
	
	def __new__(cls, x, y=None, scale=GRID_SPACING):
		"""Coord instances scale grid units (that trees use) up to SVG units (that the drawing uses). By default, the scale factor is :py:data:`~thompson.constants.GRID_SPACING` (40).
		
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
			except Exception:
				pass
			
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


def new_drawing(filename='output.svg', num_grid_lines=0):
	"""Creates a new SVG :class:`Drawing <svgwrite:svgwrite.drawing.Drawing>`. A :class:`Group <svgwrite:svgwrite.container.Group>` called *canvas* is created within the drawing, positioned slightly away from the origin. The grid that :class:`Coord` describes can optionally be drawn.
	
	CSS Selectors:
	
	- The canvas has id 'canvas'.
	- If drawn, the grid is contained in a group with class 'grid'.
		- Axis lines have class 'axis'.
		- Every fifth line has class 'major'.
	
	:param str filename: The filename associated with the drawing. **NB**: no file is created on disk until :py:meth:`dwg.save <svgwrite:svgwrite.drawing.Drawing.save>` is called.
	:param int num_grid_lines: Draws the given number of lines of the grid coordinate system, as well as axes.
	:return: A 2-tuple (drawing, canvas).
	
	.. todo::
		
		Deploy the stylesheet with the SVG
	"""
	dwg = svgwrite.Drawing(filename)
	dwg.add_stylesheet('styles.css',title=None)
	
	create_markers(dwg)
	
	canvas = svgwrite.container.Group(id='canvas')
	canvas.translate(Coord(1, 1))
	dwg.add(canvas)
	
	if num_grid_lines > 0:
		grid = svgwrite.container.Group(class_='grid')
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

def create_markers(dwg):
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

def offset_group(g, offset=Coord(0.5, 0.5)):
	"""Translates a :class:`Group <svgwrite:svgwrite.container.Group>`. If the group has a size attribute, a rectangle with class 'debug' is added, representing the size of the group."""
	g.translate(offset)
	if hasattr(g, 'size'):
		r = svgwrite.shapes.Rect(insert=-1*offset, size=Coord(g.size), class_='debug')
		g.elements.insert(0, r)
