"""The various parameters specify how to draw trees and so on. These deal with creating and positioning elemetnts; as much styling as possible is achieved via CSS."""
from svgwrite import cm, em, ex

#:When drawing trees we assign integer coordinates to each node, thus positioning the tree on a grid.
#:One grid unit is equal to this amount of `SVG units`_.
GRID_SPACING = 40

#:The radius of a tree node in `SVG units`_.
NODE_RADIUS = GRID_SPACING//2.1

#:The length of an arrow representing a group element, in grid units.
#:
#:.. todo::
#:	
#:	link to TreePair
ARROW_LENGTH = 2

#: When viewing a group element as a bijection of the interval :math:`[0, 1]`, this is the length of the interval in grid units.
GRAPH_SCALE_FACTOR = 10