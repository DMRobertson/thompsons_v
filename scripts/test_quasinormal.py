from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

qnb = example_5_3.quasinormal_basis()
p, i = example_5_3._partition_basis(qnb)
print(p)
print(example_5_3.free_factor(p, infinite=False))