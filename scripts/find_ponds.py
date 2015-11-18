from scripts import setup_script, find_examples_passing
setup_script(__file__)

from thompson.examples import random_automorphism

"""A script that hunts for ponds."""

def has_pond(aut):
	if not aut.pond_banks:
		return ''
	
	out = '{} ponds:\n\t'.format(len(aut.pond_banks))
	ponds = '\n\t'.join(  ',  '.join(str(x) for x in triple)  for triple in aut.pond_banks)
	return out + ponds

if __name__ == "__main__":
	find_examples_passing(has_pond, random_automorphism, 'find_ponds', description = "This test records examples of automorphisms which have ponds.")