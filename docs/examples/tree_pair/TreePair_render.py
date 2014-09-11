from thompson.tree_pair import TreePair
from thompson.drawing import *

x = TreePair("11000", "10100", "1 2 3")

g = x.render()
dwg, canvas = new_drawing(filename='TreePair_render.svg')
canvas.add(g)
dwg['width'], dwg['height'] = g.size + Coord(2, 2)
dwg.save()
 
g = x.render_bijection()
dwg, canvas = new_drawing(filename='TreePair_render_bijection.svg')
canvas.add(g)
dwg['width'], dwg['height'] = g.size + Coord(2, 2)
dwg.save()