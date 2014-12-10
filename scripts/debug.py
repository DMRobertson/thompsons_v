from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

import cProfile

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

from thompson import Automorphism
from thompson.examples import *

c_profile = False
filename = 'QNB_size_206'
aut = Automorphism.from_file(filename + '.aut')

if c_profile:
	print('profiling with cProfile')
	cProfile.run('X = aut.quasinormal_basis()', filename=filename+'.profile')
else:
	print('profiling with PyCallGraph')
	with PyCallGraph(output=GraphvizOutput()):
		X = aut.quasinormal_basis()

print(X)
for x in X:
	ctype, _, _ = aut.orbit_type(x, X)
	print(x, ctype)