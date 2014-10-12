import svgwrite

from thompson.constants import *
from thompson.drawing import *

g1 = svgwrite.container.Group()
g2 = svgwrite.container.Group()
coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
for centre in coords:
	circle = svgwrite.shapes.Circle(Coord(centre), r=NODE_RADIUS)
	g1.add(circle)
	g2.add(circle)
set_size(g1, Coord(2, 2),)
set_size(g2, Coord(2, 2), offset=Coord(0.5, 0.5))
g2.translate(Coord(3, 0))

dwg, canvas = new_drawing('set_size_demo.svg', num_grid_lines=5, debug_css=True)
dwg['width'], dwg['height'] = Coord(7, 7)

canvas.add(g1)
canvas.add(g2)
dwg.save()