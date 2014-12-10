from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson.automorphism import Automorphism
from thompson.orbits import *
# from thompson.examples import *
from thompson.word import Word
from pprint import pprint

import cProfile

aut = Automorphism.from_file('QNB_size_344.aut')
print('profiling the QNB')
cProfile.run('X = aut.quasinormal_basis()', filename='QNB_size_344.profile')
print(X)

for x in X:
	ctype, _, _ = aut.orbit_type(x, X)
	print(x, ctype)