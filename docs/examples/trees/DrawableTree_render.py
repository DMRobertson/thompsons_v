import random
from thompson.trees import DrawableTree, random_tree
from thompson.drawing import *

x = random_tree()
x.layout()
g = x.render()
dwg, canvas = new_drawing(filename='DrawableTree_render.svg') 
dwg['width'], dwg['height'] = g.size + Coord(2, 2)
canvas.add(g)
dwg.save()