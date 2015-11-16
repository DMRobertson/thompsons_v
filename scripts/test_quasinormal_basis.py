from scripts import setup_script
setup_script(__file__)

"""Repeatedly compute the QNBs of random automorphisms, and record the automorphism which causes an error."""

from thompson import *
from thompson.examples import *

i = 0
while True:
    i += 1
    if i % 50 == 0:
        print(i)
    aut = random_automorphism()
    try:
        aut.quasinormal_basis
    except BaseException:
        aut.save_to_file('problem.aut')
        raise