from scripts import setup_script
setup_script(__file__)

from svgwrite           import cm
from svgwrite.drawing   import Drawing
from svgwrite.shapes    import Line, Rect
from svgwrite.container import Group

default_config = {
	'scale': 400, #the width the interval [0, 1] is drawn as
}

def new_drawing(filename='plot.svg', config=default_config):
	dwg = Drawing(filename='plot.svg')
	dwg.add_stylesheet('plot.css', title='Styles for plotting elements of V_{n,r}')
	#todo eventually: embed this into the svg
	
	canvas = Group(id='canvas')
	canvas.scale(config['scale'])
	canvas.translate((0.1, 0.1))
	dwg.add(canvas)
	
	canvas.add(Rect(class_='debug'))
	return dwg, canvas

def plot(self, config=default_config):
	dwg, canvas = new_drawing(config=config)
	
	for d, r in zip(self.domain, self.range):
		x0, x1 = (float(x) for x in d.as_interval())
		y0, y1 = (float(y) for y in r.as_interval())
	
	x_axis = Line((0, 0), (1, 0), class_="axis")
	y_axis = Line((0, 0), (0, 1), class_="axis")
	canvas.add(x_axis)
	canvas.add(y_axis)
	
	dwg.save()

if __name__ == "__main__":
	from thompson.examples import standard_generator
	plot(standard_generator(0))

