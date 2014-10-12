from test import setup_script
dwg, canvas = setup_script(__file__)

"""Generates two random tree pairs *f* and *g* and forms the products ``f * g`` and ``g * f``. All four trees are output to an SVG."""

from thompson.tree_pair import *
from thompson.drawing import *

f = random_pair()
g = random_pair()
render_products(canvas, f, g)
dwg['width'], dwg['height'] = canvas.size + Coord(2, 2)
dwg.save()