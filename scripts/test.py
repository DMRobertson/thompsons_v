import sys
from pathlib import Path

"""Naughty hack to make the thompson package importable when double-clicking a python file in the /scripts directory. Adds the package root directory to sys.path.

TODO: because I'm using pathlib this implicitly relies on the pathlib package or Python 3.4. In fact maybe I should just make this a package for 3.4?
"""
def setup_script(__file__):
	"""Pass this function the __file__ argument. It will allow the script to import the thompson package and return a drawing, canvas pair to be used for drawing output."""
	path = Path(__file__)
	try:
		parent = str(path.parents[1])
	except IndexError:
		parent = str(Path('..').resolve())
	sys.path.insert(0, parent)
	# svg_file = path.stem + ".svg"
	# from thompson.drawing import new_drawing
	# return new_drawing(filename=svg_file)
