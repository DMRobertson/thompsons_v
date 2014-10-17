from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

from thompson.examples import example_4_5
basis = Generators.standard_basis(2, 1)
basis.minimal_expansion(example_4_5) == example_4_5.domain