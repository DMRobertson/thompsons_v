from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

u = Word("x a1 a1 x a1 a2 a1 x a1 a2 a2 a1 L L", 2, 1)
v1 = Word("x a1 a1 x a2 a2 x a2 a1 L L", 2, 1)
v2 = Word("x a1 a1 x a1 a2 a1 x a1 a1 L L", 2, 1)
print(cyclic_order_six.share_orbit(u, v1))
# {..., -1, 5, 11, 17, 23, ...}
print(cyclic_order_six.share_orbit(u, v2))
# {}