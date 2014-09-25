from test_script import import_thompson
dwg, canvas = setup_script(__file__)

from thompson.tree_pair import *
from thompson.drawing import *

f = random_pair()
g = random_pair()
render_products(canvas, f, g)
dwg['width'], dwg['height'] = canvas.size + Coord(2, 2)
dwg.save()