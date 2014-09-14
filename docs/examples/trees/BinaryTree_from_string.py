from thompson.trees import DrawableTree
from thompson.drawing import *

branch_patterns = "100", "1100100", "11000"
dwg, canvas = new_drawing('BinaryTree_from_string.svg')

offset = 0
height = 0
for pattern in branch_patterns:
	x = DrawableTree.from_string(pattern)
	x.layout()
	g = x.render()
	g.translate(offset)
	canvas.add(g)
	offset += g.size.x + Coord.length(1)
	height = max(height, g.size.y)

dwg['width'], dwg['height'] = Coord.unscaled(offset, height) + Coord(1, 2)
dwg.save()
