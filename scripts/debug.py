from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson.automorphism import Automorphism
from thompson.orbits import *
from thompson.examples import *
from thompson.word import Word
from pprint import pprint

psi = Automorphism.from_file('psi.aut')
phi = Automorphism.from_file('phi.aut')

rho = psi.test_conjugate_to(phi)
assert rho is not None
assert psi*rho == rho*phi