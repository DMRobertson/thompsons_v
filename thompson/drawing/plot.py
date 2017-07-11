from collections import defaultdict

from ..generators   import Generators
from ..homomorphism import Homomorphism
from ..word         import Signature

import pkg_resources

PADDING = 40
SCALE   = 400

def new_drawing(filename='plot.svg'):
	from svgwrite.drawing   import Drawing
	from svgwrite.container import Group
	size = 2 * PADDING + SCALE
	dwg = Drawing(filename, size=(size, size))
	dwg.viewbox(minx=0, miny=0, width=size, height=size)
	
	add_stylesheet(dwg)
	
	canvas = Group(id='canvas')
	canvas.translate(PADDING, PADDING)
	canvas.scale(1, -1)
	canvas.translate(0, -SCALE)
	dwg.add(canvas)
	
	return dwg, canvas

def add_stylesheet(dwg, filename='plot.css'):
	from svgwrite.container import Style
	path = pkg_resources.resource_filename('thompson.drawing', filename)
	with open(path) as f:
		contents = f.read()
	style = Style(contents)
	dwg.defs.add(style)

def include_markers(dwg, endpoints):
	from svgwrite.container import Marker
	from svgwrite.shapes    import Circle, Line, Polyline
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
	
	if not endpoints:
		return
	
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
	if signature is None:
		signature = Signature(2, 1)
	
	if levels is None:
		levels = 4
		if signature.alphabet_size > 1:
			levels -= 1
		if signature.arity > 3:
			levels -= 1
	
	from svgwrite.shapes    import Line
	done = set()
	basis = Generators.standard_basis(signature)
	lines = defaultdict(list)
	
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

def plot(*auts, filename='plot.svg', diagonal=True, endpoints=False):
	assert len({aut.signature for aut in auts if isinstance(aut, Homomorphism)}) <= 1
	try:
		signature = next(aut.signature for aut in auts if isinstance(aut, Homomorphism))
	except:
		signature = None
	
	from svgwrite.path      import Path
	from svgwrite.shapes    import Line, Polyline
	from svgwrite.container import Group
	dwg, canvas = new_drawing(filename)
	include_markers(dwg, endpoints)
	draw_grid(canvas, signature)
	
	x_axis = Polyline([(0, 0), (SCALE, 0)], class_="axis")
	y_axis = Polyline([(0, 0), (0, SCALE)], class_="axis")
	canvas.add(x_axis)
	canvas.add(y_axis)
	
	if diagonal:
		diag = Line((0,0), (SCALE, SCALE), class_="grid depth0")
		canvas.add(diag)
	
	for i, aut in enumerate(auts):
		group = Group(class_="graph", id='graph_' + str(i))
		canvas.add(group)
		last = (None, None)
		for (x0, y0, x1, y1) in graph_segments(aut):
			if last != (x0, y0):
				graph = Path(class_='graph_segment')
				group.add(graph)
				graph.push('M', SCALE * x0, SCALE * y0)
			graph.push('L', SCALE * x1, SCALE * y1)
			last = (x1, y1)
		
	dwg.save()

def graph_segments(aut):
	if isinstance(aut, Homomorphism):
		for d, r in zip(aut.domain, aut.range):
			x0, x1 = (float(x) for x in d.as_interval())
			y0, y1 = (float(y) for y in r.as_interval())
			yield (x0, y0, x1, y1)
	else:
		#this doesn't detect jumps 1 -> 0, but it's a crude starting point
		from plmaps.util import pairwise
		for (d0, r0), (d1, r1) in pairwise(aut):
			d0 = crude_mod_one(d0)
			r0 = crude_mod_one(r0)
			d1 = crude_mod_one(d1)
			r1 = crude_mod_one(r1)
			yield (d0, r0, d1, r1)

def crude_mod_one(x):
	x = float(x)
	if x > 1:
		x %= 1
	return x