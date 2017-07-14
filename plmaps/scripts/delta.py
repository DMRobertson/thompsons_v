from thompson import *
from plmaps   import *
from plmaps.debug import debug; debug()
from plmaps.util  import dump

#1. Load the definition of delta
gamma = CPL2(
	[0, F(1,4), F(3,8), F(1,2), F(9,16) , F(11,16), F(3,4)  , F(7,8), 1],
	[0, F(1,8), F(1,4), F(1,2), F(17,32), F(9,16) , F(11,16), F(3,4), 1]
)
rotation = CPL2.rotation(F(1, 2))
delta = gamma * rotation
dump(delta, "delta")

print("delta commutes with rotation by 1/2:", delta.commutes(rotation))

#2. Compute its qth power, restricted to a fundamental domain.
delta_squared = delta * delta
dump(delta_squared, "delta_squared")
dsr = delta_squared.restriction(0, F(1, 2))
dump(dsr, "dsr")

#2A. Check what DSR actually is (cf gamma.py)
x0 = PL2.from_aut(standard_generator())
y  = PL2(
	[F(0), F(1, 8), F(3,8), F(1,2), F(3,4), F(1,1)],
	[F(0), F(1,16), F(1,8), F(3,8), F(1,2), F(1,1)]
)
rescaler1 = PL2(
	[0, 1],
	[0, F(1,2)]
)
x0 ^= rescaler1
y  ^= rescaler1
assert dsr == x0 * y

#3. dsr is a one bump function, so has a cyclic centraliser.
#   construct a generator.
cent_generator = dsr.one_bump_cent_gen()
dump(cent_generator, "cent_generator")
assert cent_generator == dsr
