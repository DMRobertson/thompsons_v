from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson.examples.random import random_infinite_automorphism

while True:
	phi = random_infinite_automorphism()
	print(phi)
	