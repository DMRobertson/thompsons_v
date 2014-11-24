from test import setup_script
setup_script(__file__)

from thompson.examples import *
from thompson.automorphism import *
from pprint import pprint


psi = Automorphism.from_file('first_pond_example_psi.aut')
print(psi)
X = psi.quasinormal_basis()
print(*psi.semi_infinite_end_points())
print(X)