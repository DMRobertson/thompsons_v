from test import setup_script
setup_script(__file__)

from thompson.word import Word
from thompson.generators import Generators
from thompson.automorphism import Automorphism


leaves = ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2 a1", "x a1 a2 a2 a2", "x a2 a1", "x a2 a2"]
domain = Generators(2, 1, leaves)
range = Generators(2, 1, [leaves[0], leaves[3], leaves[4], leaves[5], leaves[2], leaves[1]])
Automorphism.reduce(domain, range)
for d, r in zip(domain, range):
	print(d, '->', r)

"""
Y = Generators(2, 1, ["x a1",       "x a2 a1",    "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
W = Generators(2, 1, ["x a1 a1 a1", "x a1 a1 a2", "x a1 a2",       "x a2 a1",       "x a2 a2"   ])

phi = Automorphism(2, 1, Y, W)
print(phi._minimal_expansion())
"""
