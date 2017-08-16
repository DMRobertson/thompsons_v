from thompson import *
from plmaps   import *
from plmaps.debug import debug; debug()
from plmaps.util  import dump, pairwise

#1. Load the definition of delta
x0 = PL2.from_aut(standard_generator(0))
x1 = PL2.from_aut(standard_generator(1))

def Sigma(s):
	return PL2( [0, 1], Word(s, (2,1)).as_interval() )

Sigma0  = Sigma('0')
Sigma00 = Sigma('00')
Sigma01 = Sigma('01')
Sigma1  = Sigma('1')
Sigma10 = Sigma('10')
Sigma11 = Sigma('11')
rot = CPL2(
	[0, 1],
	[F(1,2), F(3,2)]
)

assert (rot * rot).is_identity()

def glue(*auts, cls=PL2):
	D = list( auts[0].domain )
	R = list( auts[0].range  )
	for prev, next in pairwise(auts):
		assert next.domain[0] == prev.domain[-1]
		assert next.range [0] == prev.range [-1]
		D += next.domain[1:]
		R += next.range[1:]
	return cls(D,R)

epsilon = glue(x1^Sigma0, x1^Sigma1, cls=CPL2)
zeta = epsilon * rot
assert zeta.commutes(rot)

zeta_squared = zeta * zeta
zsr = zeta_squared.restriction(0, F(1, 2))

generators = [
	glue( x0^Sigma00, PL2.identity(F(1,4), F(1,2)) ),
	glue( x1^Sigma00, PL2.identity(F(1,4), F(1,2)) ),
	x1^Sigma0,
]
#extend to elements of C_T(zeta)_0
z0 = zeta.restriction(0, F(1,2))

for i, _ in enumerate(generators):
	generators[i] = glue( generators[i], generators[i]^z0, cls=CPL2 )
	assert generators[i].commutes(zeta)

