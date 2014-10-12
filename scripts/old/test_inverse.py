from test import setup_script
dwg, canvas = setup_script(__file__)

"""Generates many random tree pairs *f* and ensures that ``f * ~f`` is the identity."""

from thompson.tree_pair import *
from thompson.drawing import *

id = TreePair("0", "0")

def test(f, g):
	return f * g == g * f == id

def dump(f, g):
	render_products(canvas, f, ~f)
	dwg['width'], dwg['height'] = canvas.size + Coord(2, 2)
	dwg.save()

for i in range(10**4):
	if i % 100 == 0:
		print(i)
	f = random_pair()
	g = ~f
	if not test(f, g):
		dump(f, g)
		raise ValueError("Inversion failed for {!r}".format(f))

input('Test complete - no errors found. Press enter to exit.')