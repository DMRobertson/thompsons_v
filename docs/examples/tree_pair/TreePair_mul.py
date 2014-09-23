import svgwrite
from thompson.tree_pair import TreePair
from thompson.drawing import *

def layout(canvas, a, b, c, d, names=('a', 'b', 'c', 'd')):
	a = a.render(name = names[0])
	b = b.render(name = names[1])
	c = c.render(name = names[2])
	d = d.render(name = names[3])
	
	x = Coord.unscaled(max(a.size.x, c.size.x), 0) + Coord(1, 0)
	y = Coord.unscaled(0, max(a.size.y, b.size.y)) + Coord(0, 1)
	
	b.translate(x)
	c.translate(y)
	d.translate(x + y)
	
	canvas.add(a)
	canvas.add(b)
	canvas.add(c)
	canvas.add(d)
	
	return d.size + x + y + Coord(2, 2)
	

f = TreePair("100", "100", "2 1")
g = TreePair("10100", "10100", "1 3 2")
dwg, canvas = new_drawing(filename='TreePair_mul.svg')

size = layout(canvas, f, g, f*g, g*f, ['f', 'g', 'fg', 'gf'])
dwg['width'], dwg['height'] = size
dwg.save()