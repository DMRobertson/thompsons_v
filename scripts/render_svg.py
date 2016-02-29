from scripts import setup_script
setup_script(__file__)

import sys
from collections import defaultdict

from thompson import Generators

from svgwrite.drawing   import Drawing
from svgwrite.path      import Path
from svgwrite.shapes    import Circle, Line, Polygon, Polyline, Rect
from svgwrite.container import Group, Marker

PADDING = 40
SCALE   = 400

def new_drawing(filename='plot.svg'):
	dwg = Drawing(filename)
	size = 2 * PADDING + SCALE
	dwg.viewbox(minx=0, miny=0, width=size, height=size)
	
	dwg.add_stylesheet('plot.css', title='Styles for plotting elements of V_{n,r}')
	#todo eventually: embed this into the svg
	
	canvas = Group(id='canvas')
	canvas.translate(PADDING, PADDING)
	canvas.scale(1, -1)
	canvas.translate(0, -SCALE)
	dwg.add(canvas)
	
	include_markers(dwg)
	return dwg, canvas

def include_markers(dwg):
	arrowhead = Marker(
		insert=(4, 4), size=(4, 8), orient='auto', markerUnits='strokeWidth',
		id='arrowhead'
	)
	arrowhead.add(Polyline([(0, 8), (4, 4), (0, 0)]))
	dwg.defs.add(arrowhead)
	
	#TICKMARK
	tickmark = Marker(
		insert=(0.5, 5), size=(1, 10), orient='auto', markerUnits='strokeWidth',
		id='tickmark'
	)
	tickmark.add(Line((0, 0), (0, 10)))
	dwg.defs.add(tickmark)
	
	#DISCONTINUITY
	discontinuity = Marker(
		insert=(3.5, 3.5), size=(7, 7), markerUnits='strokeWidth',
		id='discontinuity'
	)
	discontinuity.add(Circle((3.5, 3.5), 2))
	dwg.defs.add(discontinuity)
	
	#CONTINUITY
	continuity = Marker(
		insert=(3.5, 3.5), size=(7, 7), markerUnits='strokeWidth',
		id='continuity'
	)
	continuity.add(Circle((3.5, 3.5), 2))
	dwg.defs.add(continuity)

def draw_grid(canvas, signature, levels=None):
	done = set()
	basis = Generators.standard_basis(signature)
	lines = defaultdict(list)
	
	if levels is None:
		levels = 4
		if signature.alphabet_size > 1:
			levels -= 1
		if signature.arity > 3:
			levels -= 1
	
	for d in range(levels):		
		for word in basis:
			end = word.as_interval()[1]
			if end in done:
				continue
			done.add(end)
			end = float(end)
			vert = Line( (SCALE * end, SCALE * 0), (SCALE * end, SCALE * 1),
				class_='grid depth' + str(d))
			hor  = Line( (SCALE * 0, SCALE * end), (SCALE * 1, SCALE * end),
				class_='grid depth' + str(d))
			lines[d].append(vert)
			lines[d].append(hor)
		
		for i in reversed(range(len(basis))):
			basis.expand(i)
		
	for d in reversed(sorted(lines)):
		for line in lines[d]:
			canvas.add(line)

def plot(self):
	dwg, canvas = new_drawing()
	draw_grid(canvas, self.signature)
	
	x_points = [(0, 0)]
	y_points = [(0, 0)]
	graph_segments = []
	for d, r in zip(self.domain, self.range):
		x0, x1 = (float(x) for x in d.as_interval())
		y0, y1 = (float(y) for y in r.as_interval())
		
		x_points.append( (SCALE * x1, 0) )
		y_points.append( (0, SCALE * y1) )
		graph_segments.append( (x0, y0, x1, y1) )
	
	x_points = sorted(x_points)
	y_points = sorted(y_points)
	
	x_axis = Polyline(x_points, class_="axis")
	y_axis = Polyline(y_points, class_="axis")
	canvas.add(x_axis)
	canvas.add(y_axis)
	
	last = (None, None)
	for (x0, y0, x1, y1) in graph_segments:
		if last != (x0, y0):
			graph = Path(class_='graph')
			canvas.add(graph)
			graph.push('M', SCALE * x0, SCALE * y0)
		graph.push('L', SCALE * x1, SCALE * y1)
		last = (x1, y1)	
	dwg.save()

if __name__ == "__main__":
	from thompson.examples import random_automorphism
	from random import randint
	x = random_automorphism(num_expansions=randint(1, 5))
	plot(x)
	print(x.signature)

