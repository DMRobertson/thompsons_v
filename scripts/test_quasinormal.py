from test import setup_script
setup_script(__file__)

from thompson.word import Word, from_string, standardise, format
from thompson.generators import Generators
from thompson.automorphism import Automorphism
# from thompson.examples import *

print(Word('x1 a1 x1 a2 x1 a3 L a2 a1 a3 x1 a2 a1 x2 a1 x1 a1 a1 x1 a1 a2 x2 L x1 a2 L a2 a1 a1 L a3 a3 a2', 3, 2))

# w = Word('x1 a1 a2', 2, 1)
# print(example_4_25[w])