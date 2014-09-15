import random
from thompson.trees import DrawableTree
from thompson.drawing import *

"""Build a string pattern of branches and nodes (ones and zeroes) by randomly joining subtrees together.
I first saw this implemented in Roman Kogan's nvTrees applet (NvTree.java).
http://www.math.tamu.edu/~romwell/nvTrees/
"""
NUM_NODES = random.randint(5, 15)
subtrees = ["0"]*NUM_NODES

N = len(subtrees)
indices = list(range(N))

for i in range(N-1):
	i, j = random.sample(indices, 2)
	subtrees[i] = "1" + subtrees[i] + subtrees[j]
	del indices[-1], subtrees[j]
	N -= 1

x = DrawableTree.from_string(subtrees[0])
x.layout()
g = x.render()
dwg, canvas = new_drawing(filename='DrawableTree_render.svg') 
dwg['width'], dwg['height'] = g.size+ Coord(2, 2)
canvas.add(g)
dwg.save()