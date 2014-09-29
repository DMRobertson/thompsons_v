from test import setup_script
setup_script(__file__)

from thompson.word import Word
from thompson.generators import Generators
from thompson.automorphism import Automorphism

Y = Generators(2, 1, ["x a1",       "x a2 a1",    "x a2 a2 a1 a1", "x a2 a2 a1 a2", "x a2 a2 a2"])
W = Generators(2, 1, ["x a1 a1 a1", "x a1 a1 a2", "x a1 a2",       "x a2 a1",       "x a2 a2"   ])

phi = Automorphism(2, 1, Y, W)

#1. Reduce the automorphism (eliminate carets)
#2. Find elements above (Y union W)
#3. Expand std basis until it contains all the elements from step 2.
#4. Test each element of this expanded basis to see if the have nice orbits.
	#a. If an element doesn't, expand it and try again with its children.
#5. When we are done the result is the smallest possible semi-normal basis for phi, i.e. a quasi-normal
