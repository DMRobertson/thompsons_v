from plmaps import *
from fractions import Fraction
psi = PL2(
	[Fraction(1, 6), Fraction(1, 2), Fraction(5, 8), Fraction(2, 3)],
	[Fraction(1, 6), Fraction(1, 4), Fraction(1, 2), Fraction(2, 3)]
)
from plmaps.debug import debug; debug()

x = psi.one_bump_cent_gen()
assert x == psi