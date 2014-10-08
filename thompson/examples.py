"""A list of examples for use in doctests."""

from .automorphism import Automorphism
from .generators import Generators

__all__ = ["cyclic_order_six", "arity_three_order_inf",
           "example_4_5", "example_4_11", "example_4_12", "example_4_25"]

#Example 4.5
domain = Generators(2, 1, ["x1 a1 a1 a1", "x1 a1 a1 a2", "x1 a1 a2",    "x1 a2 a1", "x1 a2 a2"])
range  = Generators(2, 1, ["x1 a1 a1",    "x1 a1 a2 a1", "x1 a1 a2 a2", "x1 a2 a2", "x1 a2 a1"])
example_4_5 = Automorphism(2, 1, domain, range)
#See example e.18 for characteristics

#Example 4.11
#TODO. I think this is a special element of one of Thompson's groups.
#I should put those examples in the automorphism module.
domain = Generators(2, 1, ["x a1 a1", "x a1 a2", "x a2"])
range  = Generators(2, 1, ["x a1",    "x a2 a1", "x a2 a2"])
example_4_11 = Automorphism(2, 1, domain, range)

#Example 4.12
domain = Generators(2, 1, ["x a1 a1", "x a1 a2", "x a2"])
range  = Generators(2, 1, ["x a2 a2", "x a2 a1", "x a1"])
example_4_12 = Automorphism(2, 1, domain, range)

#Example 4.25
#Todo. Characteristics?
domain = Generators(2, 1, ["x a1",       "x a2 a1",    "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
range  = Generators(2, 1, ["x a1 a1 a1", "x a1 a1 a2", "x a2 a2",       "x a1 a2",       "x a2 a1"   ])
example_4_25 = Automorphism(2, 1, domain, range)

#A completely periodic element of order six
leaves = ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2 a1", "x a1 a2 a2 a2", "x a2 a1", "x a2 a2"]
domain = Generators(2, 1, leaves)
range = Generators(2, 1, [leaves[0], leaves[3], leaves[4], leaves[5], leaves[2], leaves[1]])
cyclic_order_six = Automorphism(2, 1, domain, range)

#An example with arity n=3. 
d = Generators.standard_basis(3, 1).expand(0).expand(2).expand(0)
r = Generators.standard_basis(3, 1).expand(0).expand(2).expand(4)
r[-3:] = reversed(r[-3:])
arity_three_order_inf = Automorphism(3, 1, d, r)