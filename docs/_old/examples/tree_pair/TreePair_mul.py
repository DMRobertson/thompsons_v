import svgwrite
from thompson.tree_pair import TreePair, render_products
from thompson.drawing import *

f = TreePair("100", "100", "2 1")
g = TreePair("10100", "10100", "1 3 2")
dwg, canvas = new_drawing(filename='TreePair_mul.svg')

render_products(canvas, f, g)
dwg['width'], dwg['height'] = canvas.size
dwg.save()