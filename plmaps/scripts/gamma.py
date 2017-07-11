"""
I'm going to define an element of F of the form
    |  y 
---------
 x0 |    
where y is conjugate to x0 in F.
"""

from thompson import *
from thompson.drawing import plot
from plmaps   import *
from plmaps.util import dump
from plmaps.debug import debug; debug()

x0 = PL2.from_aut(standard_generator())
y  = PL2(
	[F(0), F(1, 8), F(3,8), F(1,2), F(3,4), F(1,1)],
	[F(0), F(1,16), F(1,8), F(3,8), F(1,2), F(1,1)]
)

result = x0.one_bump_test_conjugate(y)

rescaler1 = PL2(
	[0, 1],
	[0, F(1,2)]
)
x0 ^= rescaler1

rescaler2 = PL2(
	[0, 1],
	[F(1,2), 1]
)

y ^= rescaler2

gamma = glue(x0, y)
dump(gamma, "gamma")

#Next: construct omega and check that it commutes with gamma
phi = PL2(
	[0, F(1,4), F(1,2), F(3,4), 1],
	[0, F(1,8), F(3,8), F(1,2), 1]
)
rot = CPL2.rotation(F(1,2))
omega = glue(phi^rescaler1, ~phi^rescaler2, cls=CPL2) * rot
dump(omega, "omega")

print("gamma commutes with omega:", gamma.commutes(omega))
dump(omega, "omega")
print("gamma commutes with rotby 1/2:", gamma.commutes(rot))

print("omega^2 = id:", (omega * omega).is_identity())

#generators of the centraliser
gen1 = PL2(
	[x for x in gamma.domain if x <= F(1,2)] + [1],
	[x for x in gamma.range  if x <= F(1,2)] + [1],
)
gen2 = PL2(
	[0] + [x for x in gamma.domain if x >= F(1,2)],
	[0] + [x for x in gamma.range  if x >= F(1,2)],
)

print("Centraliser generators:")
dump(gen1, "gen1")
dump(gen2, "gen2")

print("omega commutes w/ gamma:", omega.commutes(gamma))
print("gen1  commutes w/ gamma:", omega.commutes(gamma))
print("gen2  commutes w/ gamma:", omega.commutes(gamma))
print("gen1  commutes w/ gen2 :", gen1 .commutes(gen2 ))
print("Centraliser relations")
print("gen1^omega == gen2:", gen1^omega == gen2)
print("gen2^omega == gen1:", gen2^omega == gen1)