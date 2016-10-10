from collections import defaultdict

from .automorphism import Automorphism
from .generators   import Generators

class SCSystem:
	def __init__(self, list1, list2):
		if not all( _aut_in_V(x) for x in list1 ):
			raise TypeError("An element of list1 is not an element of V(2,1)")
		if not all( _aut_in_V(x) for x in list2 ):
			raise TypeError("An element of list2 is not an element of V(2,1)")
		if len(list1) != len(list2):
			raise ValueError("List lengths differ ({} and {})".format(len(list1), len(list2)))
		if len(list1) == 0:
			raise ValueError("Element lists are empty")
		self.list1 = list1
		self.list2 = list2
		self.length = len(list1)
		
	def __str__(self):
		return "<System of {} conjugacy equations>".format(self.length)
		
	def __iter__(self):
		yield from zip(self.list1, self.list2)
	
	def check_valid_conjugacies(self):
		for index, (aut1, aut2) in enumerate(self, start=1):
			if not aut1.is_conjugate_to(aut2):
				raise ValueError("Pair number {} is not conjguate".format(i))
	
	def align_periodic_sets(self):
		#1. Assert that period M orbits for aut1 have to go to period M orbits for aut 2
		rules = []
		for index, (aut1, aut2) in enumerate(self):
			for period in aut1.periodic_orbits:
				rule = Generators((2,1)), Generators((2,1))
				for orbit in aut1.periodic_orbits[period]:
					rule[0].extend(orbit)
				for orbit in aut2.periodic_orbits[period]:
					rule[1].extend(orbit)
				rule[0].contract()
				rule[1].contract()
				rules.append(rule)
					 
		return rules
		

def _aut_in_V(x):
	return isinstance(x, Automorphism) and x.signature == (2, 1)