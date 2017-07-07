from thompson import *
from plmaps   import *

"""In this example we take the element whose graph looks like
         x0 |    
 beta:  ---------
            | x0 
"""

def dump(aut, name):
	if isinstance(aut, Automorphism):
		print(name, aut.tikz_path(), sep="\n")
	else:
		print(name, aut, sep="\n")

#1. Load the definition of beta.
beta = CPL2(
	[  0   , F(1,4), F(3,8), F(1,2), F(3,4), F( 7,8),   1   ],
	[F(1,2), F(5,8), F(3,4), F(1,1), F(9,8), F(10,8), F(3,2)]
)
dump(beta, "beta")

#2. Compute its qth power, restricted to a fundamental domain.
beta_squared = beta * beta
dump(beta_squared, "beta_squared")
bsr = beta_squared.restriction(0, F(1, 2))
dump(bsr, "bsr")

#3. bsr is a one bump function, so has a cyclic centraliser.
#   construct a generator.
cent_generator = bsr.one_bump_cent_gen()
print("cent_generator")
print(cent_generator)
assert cent_generator * cent_generator == bsr

#4. Extend cent_generator to a map of the whole circle commuting with beta.
LHS = beta.restriction_of_range(*ends(cent_generator.domain))
prod = LHS * cent_generator * ~beta
assert prod.commutes(beta_squared.restriction(F(1, 2), 1))
#glue conj and prod together
mapping = { d: r for d, r in cent_generator }
mapping.update( pair for pair in prod )
mapping = sorted(mapping.items())
domain, range = zip(*mapping);
cent_extended = CPL2(domain, range)

print("cent_extended")
print(cent_extended)

#5. Check that cent_exnteded commutes with beta.
#   If so, the <= in the centraliser descriptions are equal.
assert cent_extended.commutes(beta)


#6. Check some computations with thompson
print("Checking some computations with `thompson'")

alpha = Automorphism.from_string("""6
(2, 1) -> (2, 1)
x1 a1 a1    -> x1 a1 a1 a1
x1 a1 a2 a1 -> x1 a1 a1 a2
x1 a1 a2 a2 -> x1 a1 a2   
x1 a2 a1    -> x1 a2 a1 a1
x1 a2 a2 a1 -> x1 a2 a1 a2
x1 a2 a2 a2 -> x1 a2 a2   """)
rot = Automorphism.rotation(1, 2)
beta2 = alpha * rot
dump(beta2, "beta2")

dump(beta2 * beta2, "beta_squared")