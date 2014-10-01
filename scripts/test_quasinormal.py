from test import setup_script
setup_script(__file__)

from thompson.word import *
from thompson.generators import Generators
from thompson.automorphism import Automorphism
from thompson.examples import *

print(lambda_arguments(from_string('x a1 a2 x a2 a2 L x a1 a1 L'), 2))

w = Word('x1 a1 a2', 2, 1)
print(example_4_25[w])