from test import setup_script
setup_script(__file__)
from pprint import pprint

from thompson.examples import random_automorphism

r"""Not sure what this was for. I think this was testing a lemma of mine that if w \in QNB is in a semi-inf component, then all words between w and the component endpoint (inclusively) belong to the QNB."""

i = 0
finished = False
while not finished:
	i += 1
	if i % 100 == 0:
		print(i)
	aut = random_automorphism()
	basis = aut.quasinormal_basis
	for gen in basis:
		ctype, images, _ = aut.orbit_type(gen, basis)
		if ctype.is_type_B():
			power = ctype.characteristic[0]
			
			if power > 0:
				indices = range(-power, 0)
			else:
				indices = range(1, abs(power)+1)
			
			for j in indices:
				try:
					x = images[j]
				except KeyError:
					continue
				if x not in basis:
					print(aut)
					print(basis)
					print(gen)
					print(ctype)
					pprint(images)
					finished = True
					break
		if finished:
			break