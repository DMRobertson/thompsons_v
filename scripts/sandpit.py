from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson import *
from thompson.examples import *

#Two sides of a pond
u = Word('x a1 a1 a1 a1 a1 a1 a1 a2', (2, 1))
v = Word('x a2 a2 a1 a1 a2', (2, 1))

print(first_pond_example_phi.quasinormal_basis())
print(first_pond_example_phi.pond_banks)

print(first_pond_example_phi.share_orbit(u, v))
