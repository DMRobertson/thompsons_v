from collections import defaultdict, namedtuple
from functools   import total_ordering

from .automorphism import Automorphism
from .generators   import Generators
from .cantorsubset import CantorSubset

__all__ = ["SCSystem", "print_rules"]

@total_ordering
class Rule:
	def __init__(self, source, target, id=None):
		if not isinstance(source, CantorSubset):
			raise TypeError("source is not a CantorSubset")
		if not isinstance(target, CantorSubset):
			raise TypeError("target is not a CantorSubset")
		self.source = source
		self.target = target
		self.id	 = id
	
	@classmethod
	def empty(cls, id=None):
		return cls(CantorSubset((2,1)), CantorSubset((2,1)), id)
	
	@classmethod
	def full(cls, id=None):
		return cls(CantorSubset.standard_basis((2,1)), CantorSubset.standard_basis((2,1)), id)
	
	def __str__(self):
		return str(self.source) + ' -> ' + str(self.target)
	
	def __iter__(self):
		yield self.source
		yield self.target
	
	def __getitem__(self, key):
		if key == 0:
			return self.source
		elif key == 1:
			return self.target
	
	def deduce(self, d, e):
		result = type(self)(d.image_of_set(self.source), e.image_of_set(self.target))
		result.source.simplify()
		result.target.simplify()
		return result
	
	def __eq__(self, other):
		return self.source == other.source and self.target == other.target
	
	def __lt__(self, other):
		return self.source < other.source or (self.source == other.source and self.target < other.target)
	
	def __and__(self, other):
		return type(self)(self.source & other.source, self.target & other.target)
	
	def __invert__(self):
		return type(self)(~self.source, ~self.target, self.id)
	
	def __bool__(self):
		return bool(self.source) or bool(self.target)

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
				raise ValueError("Pair number {} is not conjguate".format(index))
	
	def periodic_constraints(self):
		#1. Assert that period M orbits for aut1 have to go to period M orbits for aut 2
		rules = []
		for index, (aut1, aut2) in enumerate(self):
			for period in aut1.periodic_orbits:
				#a. Make a rule whose source is (aut1's period m elements)
				rule = Rule.empty()
				for orbit in aut1.periodic_orbits[period]:
					rule.source.extend(orbit)
				rule.source.simplify()
				#b. If the source is the entire Cantor set, the target will be too when we compute it.
				#   Skip this: there's no point considering a rule saying C -> C.
				if rule.source.is_entire_Cantor_set():
					continue
				#   otherwise we know we will make use of this rule
				rules.append(rule)
				#c. The target for the rule is (aut2's periodic elements)
				for orbit in aut2.periodic_orbits[period]:
					rule.target.extend(orbit)
				rule.target.simplify()
		return rules
	
	def deduce(self, rules):
		for index in range(len(rules)):
			rule = rules[index]
			rule.id = None
			#d. Impose the restrictions we get from the other conjugate pairs
			for otheraut1, otheraut2 in self:
				extra_rule = rule.deduce(otheraut1, otheraut2)
				rules.append(extra_rule)
		#TODO: make a `ruleset' class which tries to efficiently prevent duplicates
		simplify_rules(rules)
		return rules
	
	def reduce_to_atoms(self, rules):
		#3. Compute the complementary rules
		comp_rules = [~rule for rule in rules]
		
		#4. Compute the atomic intersections
		atoms = []
		for atomic_rule in self.filtered_atoms(rules, comp_rules):
			#Don't remember rules that say empty <-> empty
			if atomic_rule.source:
				atoms.append(atomic_rule)
		
		#5. Check that the union of S atoms = Cantor set. The same for T atoms
		for index in {0, 1}:
			subset = CantorSubset((2,1))
			for atom in atoms:
				subset.extend(atom[index])
			subset.simplify()
			assert subset.is_entire_Cantor_set()
		return atoms
	
	def filtered_atoms(self, rules, comp_rules):
		max_depth = len(rules) - 1
		intersection_stack = {-1: Rule.full(), 0: comp_rules[0]}
		choices_stack = [False]
		
		index = 0
		depth = 0
		force_skip = False
		while choices_stack:
			if force_skip or depth == max_depth:
				if not force_skip:
					#Check that the rule is consistant: not mapping empty <-> nonempty
					if bool(rule.source) != bool(rule.target):
						raise ValueError("Inconsistency with constraint #" + str(index))
					yield intersection_stack[depth]
				force_skip = False
				
				#b. If we're a right child:
				while choices_stack and choices_stack[-1]:
					choices_stack.pop()
					del intersection_stack[depth]
					depth -= 1
				#If we're a left child of something, move to the right child
				if choices_stack:
					index += 1 << (max_depth - depth)
					choices_stack[-1] = True
					intersection_stack[depth] = intersection_stack[depth - 1] & rules[depth]
			
			#a. If we're not at a leaf, just pick the left child (intersecting with a complement)
			elif depth < max_depth:
				depth += 1
				choices_stack.append(False)
				rule = intersection_stack[depth-1] & comp_rules[depth]
				if not rule:
					force_skip = True
				else:
					intersection_stack[depth] = rule
	
	def iterate_rules(self, rules):
		rules = self.deduce(rules)
		print_rules(rules, "Deduced new rules")
		rules = self.reduce_to_atoms(rules)
		print_rules(rules, "Reduce to atoms")
		return rules
	
	def periodic_orbit_constraints(self):
		rules = self.periodic_constraints()
		print_rules(rules, "Initial rules")
		rules = self.reduce_to_atoms(rules)
		print_rules(rules, "Reduce to atoms")
		for i in range(1):
			rules = self.iterate_rules(rules)
		return rules
	
	def coarsely_conjuate(self, atoms):
		#need atomic rules here, not just rules
		domain   = Generators(self.list1[0].signature)
		codomain = Generators(self.list2[0].signature)
		for source, target in atoms:
			size = max(len(source), len(target))
			source.expand_to_size(size)
			target.expand_to_size(size)
			domain.extend(source)
			codomain.extend(target)
		conjugator = Automorphism(domain, codomain)
		self ^= conjugator
		return conjugator
	
	def __ixor__(self, conjugator):
		if not (isinstance(conjugator, Automorphism) and conjugator.signature == self.list1[0].signature):
			raise TypeError("Can only simultaneously conjugate the system by an Automorphism of the correct signature")
		for index, _ in enumerate(self.list1):
			self.list1[index] ^= conjugator
		self.conjugator *= conjugator

def _aut_in_V(x):
	return isinstance(x, Automorphism) and x.signature == (2, 1)

def simplify_rules(rules):
	rules.sort()
	remove_duplicates_sorted(rules)
	#
	
def remove_duplicates_sorted(rules):
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

def print_rules(rules, message=None):
	if message is not None:
		print(message + ":", "(count={})".format(len(rules)))
	widths = [0, 0]
	for rule in rules:
		widths[0] = max(len(str(rule.source)), widths[0])
		widths[1] = max(len(str(rule.target)), widths[1])
	
	fmt_string = "{{:{}}} -> {{:{}}}".format(*widths)
	length = max(
		( len(bin(rule.id)) for rule in rules if isinstance(rule.id, int) )
		, default=0)
	if length > 0:
		fmt_string = "{{:0{}b}}: ".format(length) + fmt_string
		for rule in rules:
			print(fmt_string.format(rule.id or -1, rule.source, rule.target))
	else:
		for rule in rules:
			print(fmt_string.format(rule.source, rule.target))
