"""
I'm going to define an element of F of the form
    |  y 
---------
 x0 |    
where y is conjugate to x0 in F.
"""

from thompson import *
from plmaps   import *

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

#DMR: why does __ixor__ not work?
x0 = x0 ^ rescaler1