from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson import *
from thompson.examples import debug_charactersistics_phi as phi
from thompson.examples import debug_charactersistics_psi as psi

assert phi.test_conjugate_to(psi) is not None

print(phi)
phi.dump_QNB()
phi.characteristics(print=True)
print(phi.semi_infinite_end_points())
print()

print(psi)
psi.dump_QNB()
psi.characteristics(print=True)
print(psi.semi_infinite_end_points())

assert phi.characteristics() == psi.characteristics()