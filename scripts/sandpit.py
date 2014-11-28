from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson.automorphism import Automorphism
from thompson.orbits import *
# from thompson.examples import *
from thompson.word import Word
from pprint import pprint

aut = Automorphism.from_file('../thompson/examples/example_5_3.aut')
# p, i = aut.free_factors()

print(aut)
print(aut.quasinormal_basis())
# print(p)
# print(i)