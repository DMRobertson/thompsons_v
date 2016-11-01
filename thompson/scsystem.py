from collections import defaultdict

from .automorphism import Automorphism
from .generators   import Generators
from .cantorsubset import CantorSubset

__all__ = ["SCSystem", "print_rules"]

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
			
		self.check_valid_conjugacies(list1, list2)
		self.list1 = list1
		self.list1_original = list1.copy()
		self.list2 = list2
		self.conjugator = Automorphism.identity((2,1))
		self.length = len(list1)
		
		# TODO: check for pairs (d_i, e_i) and (d_j, e_j) with either entry the same
		# if they are the same pair then remove one copy;
		# if they are different pairs the system is not SC
		
	def __len__(self):
		return self.length
	
	def __str__(self):
		return "<System of {} conjugacy equations>".format(self.length)
	
	def __iter__(self):
		yield from zip(self.list1, self.list2)
	
	def __bool__(self):
		return all(aut1 == aut2 for aut1, aut2 in self)
	
	def display(self):
		conjugated = self.conjugator != 1
		for i, (aut1, aut2) in enumerate(self, start=1):
			from thompson.drawing import forest
			from IPython.display import display
			name = "$d_{}".format(i)
			if conjugated:
				name += "^c"
			name += "$"
			display( forest(aut1, name=name) )
			display( forest(aut2, name="$e_{}$".format(i)) )
	
	@staticmethod
	def check_valid_conjugacies(list1, list2):
		for index, (aut1, aut2) in enumerate(zip(list1, list2)):
			if not aut1.is_conjugate_to(aut2):
				raise ValueError("Pair number {} is not conjguate".format(i))
	
	def periodic_constraints(self):
		#1. Assert that period M orbits for aut1 have to go to period M orbits for aut 2
		rules = []
		for index, (aut1, aut2) in enumerate(self):
			for period in aut1.periodic_orbits:
				rule = (CantorSubset((2,1)), CantorSubset((2,1)))
				for orbit in aut1.periodic_orbits[period]:
					rule[0].extend(orbit)
				rule[0].simplify()
				if rule[0].is_entire_Cantor_set():
					#Have checked in the initialiser that pairs are conjugate, so rule[0] == rule[1] == entire cantor set
					continue
				for orbit in aut2.periodic_orbits[period]:
					rule[1].extend(orbit)
				rule[1].simplify()
				rules.append(rule)
		#2. Impose the restrictions we get from the other conjugate pairs relations
				for index2, (otheraut1, otheraut2) in enumerate(self):
					if index2 == index:
						#would just reproduce the rule we already have
						continue
					conjugated_rule = ( otheraut1.image_of_set(rule[0]), otheraut2.image_of_set(rule[1]) )
					conjugated_rule[0].simplify()
					conjugated_rule[1].simplify()
					rules.append(conjugated_rule)
		
		remove_duplicates(rules)
		print("Found {} rules".format(len(rules)))
		#3. Compute the complementary rules
		comp_rules = [(~set1, ~set2) for set1, set2 in rules]
		
		#4. Compute the atomic intersections
		atoms = {}
		for index, atomic_rule in self.filtered_atoms(rules, comp_rules):
			#Check that the rule is consistant: not mapping empty <-> nonempty
			if bool(atomic_rule[0]) != bool(atomic_rule[1]):
				raise ValueError("Problem with constraint #" + str(index))
			#Don't remember rules that say empty <-> empty
			if atomic_rule[0]:
				atoms[index] = atomic_rule
		
		return rules, comp_rules, atoms
		
		#5. Check that the union of S atoms = Cantor set. The same for T atoms
		for index in {0, 1}:
			subset = CantorSubset((2,1))
			for atom in atoms.values():
				subset.extend(atom[index])
			subset.simplify()
			assert subset.is_entire_Cantor_set()
				
		#6. Construct a conjugator
		return rules, comp_rules, atoms
		
	def filtered_atoms(self, rules, comp_rules):
		for index in range( 2**len(rules) ):
			bitmask = bits(index, len(rules))
			rule = [ CantorSubset.standard_basis((2,1)), CantorSubset.standard_basis((2,1)) ]
			for j, bit in enumerate(bitmask):
				next_intersectand = (rules if bit else comp_rules)[j]
				rule[0] &= next_intersectand[0]
				rule[1] &= next_intersectand[1]
			yield index, rule
	
	def align_periodic_orbits(self):
		_, _, atoms = self.periodic_constraints()
		domain = Generators(self.list1[0].signature)
		range  = Generators(self.list2[0].signature)
		for source, target in atoms.values():
			size = max(len(source), len(target))
			source.expand_to_size(size)
			target.expand_to_size(size)
			domain.extend(source)
			range.extend(target)
		conjugator = Automorphism(domain, range)
		for index, _ in enumerate(self.list1):
			self.list1[index] ^= conjugator
		self.conjugator *= conjugator
		return conjugator

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

def bits(num, num_bits):
	for _ in range(num_bits):
		yield num % 2
		num >>= 1

def print_rules(rules):
	widths = [0, 0]
	if isinstance(rules, dict):
		_rules = rules.values()
	else:
		_rules = rules
	for rule in _rules:
		widths[0] = max(len(str(rule[0])), widths[0])
		widths[1] = max(len(str(rule[1])), widths[1])
	
	fmt_string = "{{:{}}} -> {{:{}}}".format(*widths) 
	if isinstance(rules, dict):
		length = len(bin(max(rules))) - 2
		fmt_string = "{{:0{}b}}: ".format(length) + fmt_string
		for index in sorted(rules):
			rule = rules[index]
			print(fmt_string.format(index, rule[0], rule[1]))
	else:
		for rule in rules:
			print(fmt_string.format(rule[0], rule[1]))
			

