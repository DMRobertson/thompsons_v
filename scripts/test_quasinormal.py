from test import setup_script
setup_script(__file__)

from thompson.word import Word
from thompson.generators import Generators
from thompson.automorphism import Automorphism
from thompson.examples import *

w = Word('x1 a1 a2', 2, 1)
print(example_4_25[w])