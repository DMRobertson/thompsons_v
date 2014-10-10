from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

from thompson.generators import Generators
basis = Generators.standard_basis(2, 1).expand(0).expand(0).expand(0)
Word('x a2 x a1 L', 2, 1).max_depth_to(basis)
