from collections import defaultdict

from .automorphism import Automorphism
from .cantorsubset import CantorSubset

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
				rule = CantorSubset((2,1)), CantorSubset((2,1))
				for orbit in aut1.periodic_orbits[period]:
					rule[0].extend(orbit)
				for orbit in aut2.periodic_orbits[period]:
					rule[1].extend(orbit)
				rule[0].simplify()
				rule[1].simplify()
				rules.append(rule)
				#2. Impose the restrictions we get from the other conjugation relations
				for index2, (otheraut1, otheraut2) in enumerate(self):
					if index2 == index:
						continue
					conjugated_rule = ( otheraut1.image_of_set(rule[0]), otheraut2.image_of_set(rule[1]) )
					conjugated_rule[0].simplify()
					conjugated_rule[1].simplify()
					rules.append(conjugated_rule)
		
		remove_duplicates(rules)
		#3. 
		#atoms = [None] * 2 ** n
		return rules
		

def _aut_in_V(x):
	return isinstance(x, Automorphism) and x.signature == (2, 1)

def remove_duplicates(rules):
	rules.sort()
	i = 0
	while i + 1 < len(rules):
		if rules[i] == rules[i+1]:
			del rules[i+1]
		else:
		 	i += 1
	