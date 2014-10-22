from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

qnb = example_4_5.quasinormal_basis()
p, i = example_4_5._partition_basis(qnb)
print(example_4_5.free_factor(p, infinite=False))