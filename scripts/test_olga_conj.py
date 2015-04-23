from scripts import setup_script
setup_script(__file__)

"""A place to play around with the package."""

from thompson import *
from thompson.examples import *

f  = load_example('olga_f')
g  = load_example('olga_g')
gg = load_example('olga_gdash')
h  = load_example('olga_h')

rho = f.test_conjugate_to(g)
print("f is conjugate to g:", rho is not None)
#I can't remember which way round functions compose or conjugation works any more.
#I think Olga's using the opposite convention to us, though.
print("Python produced a valid conjugator:", f * rho == rho * g)
print("h is a valid conjugator:", f * h == h * g)

rho = f.test_conjugate_to(gg)
print("f is conjugate to g':", rho is not None)

print('\nQNB data for f:')
f.dump_QNB()

print('\nQNB data for g:')
g.dump_QNB()

print('\nQNB data for g\':')
gg.dump_QNB()

