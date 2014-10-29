"""A list of examples for use in doctests."""

from .word import Word
from .generators import Generators
from .automorphism import Automorphism
from .factors import PeriodicFactor, InfiniteFactor

__all__ = [
	#Made-up examples
	"cyclic_order_six", "arity_three_order_inf", "arity_four", "alphabet_size_two",
	#Paper section 4
	"example_4_5", "example_4_11", "example_4_12", "example_4_25", "example_4_25_i",
	#Conjugacy
	"example_5_3",
		#Pure Periodic
		"example_5_9", "example_5_12_psi", "example_5_12_phi", "example_5_12_rho",
		"example_5_9_p", "example_5_12_psi_p", "example_5_12_phi_p",
		#Pure Infinite
		"example_5_26"]

#TODO. Have a method that reads homomorphisms in from a file. Writing the Python syntax around all this is boring.
#TODO. Allow powers in the definition of words e.g. a1^4?

#Example 4.5
domain = Generators((2, 1), ["x1 a1 a1 a1", "x1 a1 a1 a2", "x1 a1 a2",    "x1 a2 a1", "x1 a2 a2"])
range  = Generators((2, 1), ["x1 a1 a1",    "x1 a1 a2 a1", "x1 a1 a2 a2", "x1 a2 a2", "x1 a2 a1"])
example_4_5 = Automorphism(domain, range)
#See example e.18 for characteristics

#Example 4.11
#TODO. I think this is a special element of one of Thompson's groups.
domain = Generators((2, 1), ["x a1 a1", "x a1 a2", "x a2"])
range  = Generators((2, 1), ["x a1",    "x a2 a1", "x a2 a2"])
example_4_11 = Automorphism(domain, range)

#Example 4.12
domain = Generators((2, 1), ["x a1 a1", "x a1 a2", "x a2"])
range  = Generators((2, 1), ["x a2 a2", "x a2 a1", "x a1"])
example_4_12 = Automorphism(domain, range)

#Example 4.25
domain = Generators((2, 1), ["x a1",       "x a2 a1",    "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
range  = Generators((2, 1), ["x a1 a1 a1", "x a1 a1 a2", "x a2 a2",       "x a1 a2",       "x a2 a1"   ])
example_4_25 = Automorphism(domain, range)
example_4_25_i = example_4_25.free_factor(example_4_25.quasinormal_basis(), infinite=True)

#Example 5.3
domain = ["x a1 a1 a1 a1", "x a1 a1 a1 a2", "x a1 a1 a2",    "x a1 a2 a1", "x a1 a2 a2", "x a2 a1",    "x a2 a2 a1", "x a2 a2 a2"]
range  = ["x a1 a1 a1",    "x a1 a1 a2 a1", "x a1 a1 a2 a2", "x a1 a2 a2", "x a1 a2 a1", "x a2 a1 a1", "x a2 a1 a2", "x a2 a2"]
d = Generators((2, 1), domain)
r = Generators((2, 1), range)
example_5_3 = Automorphism(d, r)

#Pure Periodic
#Example 5.9
domain = ["x a1 a1 a1", "x a1 a1 a2", "x a1 a2",    "x a2 a1 a1", "x a2 a1 a2", "x a2 a2 a1", "x a2 a2 a2"]
range  = ["x a1 a1 a2", "x a1 a2",    "x a1 a1 a1", "x a2 a1 a2", "x a2 a1 a1", "x a2 a2 a2", "x a2 a2 a1"]
d = Generators((2, 1), domain)
r = Generators((2, 1), range)
example_5_9 = Automorphism(d, r)
example_5_9_p = example_5_9.free_factor(example_5_9.quasinormal_basis())

#Example 5.12
domain = ["x a1 a1 a1 a1", "x a1 a1 a1 a2", "x a1 a1 a2", "x a1 a2",    "x a2 a1", "x a2 a2"]
range  = ["x a1 a1 a1 a2", "x a1 a1 a1 a1", "x a1 a2",    "x a1 a1 a2", "x a2 a2", "x a2 a1"]
example_5_12_psi = Automorphism(Generators((2, 1), domain), Generators((2, 1), range))
example_5_12_psi_p = example_5_12_psi.free_factor(example_5_12_psi.quasinormal_basis())

domain = ["x a1", "x a2"]
range  = ["x a2", "x a1"]
example_5_12_phi = Automorphism(Generators((2, 1), domain), Generators((2, 1), range))
example_5_12_phi_p = example_5_12_phi.free_factor(example_5_12_phi.quasinormal_basis())

domain = ["x a1 a1 a1 a1", "x a1 a1 a1 a2", "x a1 a1 a2", "x a1 a2",    "x a2 a1", "x a2 a2"]
range  = ["x a1 a1 a1",    "x a2 a1 a1",    "x a1 a1 a2", "x a2 a1 a2", "x a1 a2", "x a2 a2"]
example_5_12_rho = Automorphism(Generators((2, 1), domain), Generators((2, 1), range))
# example_5_12_rho_p = example_5_12_rho.free_factor(example_5_12_rho.quasinormal_basis())

#Example 5.26
domain = ["x a1",       "x a2 a1 a1", "x a2 a1 a2", "x a2 a2"]
range  = ["x a1 a1 a1", "x a2",       "x a1 a1 a2", "x a1 a2"]
example_5_26 = Automorphism(Generators((2, 1), domain), Generators((2, 1), range))

#Made up examples
#A completely periodic element of order six
leaves = ["x a1 a1", "x a1 a2 a1", "x a1 a2 a2 a1", "x a1 a2 a2 a2", "x a2 a1", "x a2 a2"]
domain = Generators((2, 1), leaves)
range = Generators((2, 1), [leaves[0], leaves[3], leaves[4], leaves[5], leaves[2], leaves[1]])
"""Reduced forms are 
domain = Generators((2, 1), ["x a1 a1", "x a1 a2 a1",    "x a1 a2 a2", "x a2 a1",       "x a2 a2"])
range  = Generators((2, 1), ["x a1 a1", "x a1 a2 a2 a2", "x a2",       "x a1 a2 a2 a1", "x a1 a2 a1"])
but let's leave the defn unreduced to test the _reduce method."""
cyclic_order_six = Automorphism(domain, range)

#An example with arity n=3. 
d = Generators.standard_basis((3, 1)).expand(0).expand(2).expand(0)
r = Generators.standard_basis((3, 1)).expand(0).expand(2).expand(4)
r[-3:] = reversed(r[-3:])
arity_three_order_inf = Automorphism(d, r)

#An example with arity n=4. 
d = Generators.standard_basis((4, 1)).expand(0).expand(3).expand(0)
r = ["x a1 a2", "x a1 a4", "x a1 a1", "x a1 a3", "x a3 a2", "x a3 a3", "x a3 a4", "x a2", "x a3 a1", "x a4"]
r = Generators((4, 1), r)
arity_four = Automorphism(d, r)

#Two trees! Alphabet size r=2.
d = Generators.standard_basis((3, 2)).expand(1).expand(0)
r = ["x1 a1", "x1 a2 a3", "x2", "x1 a3", "x1 a2 a2", "x1 a2 a1"]
r = Generators((3, 2), r)
alphabet_size_two = Automorphism(d, r)