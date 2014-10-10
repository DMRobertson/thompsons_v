from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

# arity_three_order_inf.dump_mapping()
# basis = example_4_12._minimal_expansion()
# basis.expand(0)
# print(basis)
# dump_orbit_types(example_4_12, basis)
dump_orbit_types(example_4_5, example_4_5.domain)