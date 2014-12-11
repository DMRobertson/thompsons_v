from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson import *
from thompson.examples import *

print(nathan_pond_example)
nathan_pond_example.quasinormal_basis()
nathan_pond_example.dump_QNB()
for (ell, k, r) in nathan_pond_example.pond_banks:
	print(ell, k, r)