from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson import *
from thompson.examples import *

sig = (3, 2)
domain = Generators(sig)
range = Generators(sig)

for line in """x1 a1 -> x2 a2
x1 a2 -> x2 a3
x1 a3 -> x1
x2    -> x2 a1""".split('\n'):
	d, r = line.split('->')
	d = Word(d.strip(), sig)
	r = Word(r.strip(), sig)
	domain.append(d)
	range.append(r)

aut = Automorphism(domain, range)
print(aut)
print(aut.cycles_order())

i = 0
while True:
	i += 1
	if i % 50 == 0:
		print(i)
	aut = random_automorphism_in_T()
	assert aut.cycles_order(), str(aut)