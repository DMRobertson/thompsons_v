import sys
from pathlib import Path

#Naughty hack to make the thompson package importable. Adds the root directory to sys.path
def import_thompson(__file__):
	"""Pass this function the __file__ argument. It will allow the script to import the thompson package and return a drawing, canvas pair to be used for drawing output."""
	path = Path(__file__)
	sys.path.insert(0, str(path.parents[1]))
	svg_file = path.stem + ".svg"
	from thompson.drawing import new_drawing
	return new_drawing(filename=svg_file)
