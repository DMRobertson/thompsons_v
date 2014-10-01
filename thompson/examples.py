"""A list of examples for use in doctests."""

from .automorphism import Automorphism
from .generators import Generators

__all__ = ["cyclic_order_six", "example_4_25"]

leaves = ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2 a1", "x a1 a2 a2 a2", "x a2 a1", "x a2 a2"]
domain = Generators(2, 1, leaves)
range = Generators(2, 1, [leaves[0], leaves[3], leaves[4], leaves[5], leaves[2], leaves[1]])
cyclic_order_six = Automorphism(2, 1, domain, range)

Y = Generators(2, 1, ["x a1",       "x a2 a1",    "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
W = Generators(2, 1, ["x a1 a1 a1", "x a1 a1 a2", "x a2 a2",       "x a1 a2",       "x a2 a1"   ])
example_4_25 = Automorphism(2, 1, Y, W)