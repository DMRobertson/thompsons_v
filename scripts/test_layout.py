"""import sys, os
sys.path.insert(0, os.path.abspath('.'))"""

from thompson.trees import DrawableTree
from thompson.drawing import new_drawing, Coord

left = DrawableTree()
left.add_child().add_child().parent.add_child(True)

right = DrawableTree.from_string("10100")

root = DrawableTree()
root.set_child(left)
root.set_child(right, right_child=True)

dwg, canvas = new_drawing(num_grid_lines=20)
def draw(node, offset=(0, 0)):
	print(node.is_strict())
	node.layout()
	g = node.render()
	g.translate(Coord(offset))
	canvas.add(g)

draw(left)
draw(right, (0, 5))
draw(root, (0, 10))


dwg.save()