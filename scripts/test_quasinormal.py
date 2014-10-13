from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

dump_orbit_types(example_4_5, example_4_5.domain)

# u  = Word('x a2 a2 a1 a1 a2', 2, 1)
# v1 = Word('x a1 a2', 2, 1)
# v2 = Word('x a1 a1 a2', 2, 1)
# print('****SHARE?', repr(example_4_25.share_orbit(u, v1)))
# None
# print('****SHARE?', repr(example_4_25.share_orbit(u, v2)))
# 3