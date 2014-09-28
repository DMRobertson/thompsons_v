from test import setup_script
setup_script(__file__)

from thompson.word import *
x = "x a1 a2 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L"
print(x)
print(Word(x, 2, 1))
# from thompson.generators import Generators
# from thompson.automorphism import Automorphism

# Y = Generators(["x a1",       "x a2 a1",    "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
# W = Generators(["x a1 a1 a1", "x a1 a1 a2", "x a1 a2",       "x a2 a1",       "x a2 a2"   ])

# phi = Automorphism(Y, W)