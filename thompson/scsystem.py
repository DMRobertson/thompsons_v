from bisect           import bisect_left
from collections      import defaultdict, namedtuple
from functools        import total_ordering
from os.path          import join, exists
from warnings         import warn
from random           import randint

from .automorphism    import Automorphism
from .generators      import Generators
from .cantorsubset    import CantorSubset
from .examples.random import random_automorphism

__all__ = ["SCSystem", "RuleSet"]

@total_ordering
class Rule:
	#todo: this should really be immutable
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
	
	def is_trivial(self):
		s_status = self.source.status()
		t_status = self.target.status()
		if s_status != t_status:
			#TODO: use our own error class here so we can catch it and only it in iterate() more easily
			raise ValueError("Incompatible source and target: {}, {}".format(self.source, self.target))
		return s_status in {-1, 1}
	
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
	
	def simplify(self):
		self.source.simplify()
		self.target.simplify()
	
	def deduce(self, d, e):
		result = type(self)(d.image_of_set(self.source), e.image_of_set(self.target))
		result.source.simplify()
		result.target.simplify()
		return result

class RuleSet:
	def __init__(self, rules=tuple()):
		self.rules = [rule for rule in rules if not rule.is_trivial()]
		self.rules.sort()
		self.atomic = False
	
	def __len__(self):
		return len(self.rules)
	
	def __iter__(self):
		return iter(self.rules)
	
	def __str__(self):
		output = str(len(self))
		if self.atomic:
			output += " atomic"
		output += " rules:\n"
		
		widths = [0, 0]
		id_length = 0
		for rule in self:
			widths[0] = max(len(str(rule.source)), widths[0])
			widths[1] = max(len(str(rule.target)), widths[1])
			if rule.id is not None:
				id_length = max( len(str(rule.id)), id_length)
		
		fmt_string = "{{:{}}} -> {{:{}}}\n".format(*widths)
		if id_length > 0:
			fmt_string = "{{:{}}}: ".format(id_length) + fmt_string
			for rule in self:
				output += fmt_string.format(rule.id or "", rule.source, rule.target)
		else:
			for rule in self:
				output += fmt_string.format(rule.source, rule.target)
		return output[:-1]
	
	def add(self, rule):
		#Return True if rule added, else False
		if rule.is_trivial():
			return False
		if len(self) == 0:
			self.rules.append(rule)
		#Is this element already present?
		index = bisect_left(self.rules, rule)
		if index >= len(self) or self.rules[index] != rule:
			#insertion is O(n); append then sort is O(log n)
			self.rules.append(rule)
			self.rules.sort()
			self.atomic = False
			return True
		return False
	
	def atomise(self):
		if self.atomic:
			return
		atoms = []
		for atomic_rule in self.filtered_atoms():
			if not atomic_rule.is_trivial():
				atoms.append(atomic_rule)
		atoms.sort()
		
		#5. Check that the union of S atoms = Cantor set. The same for T atoms
		for index in {0, 1}:
			subset = CantorSubset((2,1))
			for atom in atoms:
				subset.extend(atom[index])
			subset.simplify()
			assert subset.is_entire_Cantor_set()
		self.rules = atoms
		self.atomic = True
	
	def filtered_atoms(self, label=True):
		comp_rules = [~rule for rule in self]
		max_depth = len(self) - 1
		intersection_stack = {-1: Rule.full(), 0: comp_rules[0]}
		choices_stack = [False]
		
		fmt_string = "{{:0{}b}}".format(len(self))
		index = 0
		depth = 0
		force_skip = False
		while choices_stack:
			if force_skip or depth == max_depth:
				if not force_skip:
					#Check that the rule is consistant: not mapping empty <-> nonempty
					intersection_stack[depth].id = fmt_string.format(index)
					yield intersection_stack[depth]
				force_skip = False
				
				#b. If we're a right child:
				index += 1 << (max_depth - depth)
				while choices_stack and choices_stack[-1]:
					choices_stack.pop()
					del intersection_stack[depth]
					depth -= 1
				#If we're a left child of something, move to the right child
				if choices_stack:
					choices_stack[-1] = True
					intersection_stack[depth] = intersection_stack[depth - 1] & self.rules[depth]
			
			#a. If we're not at a leaf, just pick the left child (intersecting with a complement)
			elif depth < max_depth:
				depth += 1
				choices_stack.append(False)
				rule = intersection_stack[depth-1] & comp_rules[depth]
				if not rule:
					force_skip = True
				else:
					intersection_stack[depth] = rule
	
	def clear_labels(self):
		for rule in self.rules:
			rule.id = None
	
	def deduce_new(self, system):
		self.clear_labels()
		starting_rules = self.rules.copy()
		for rule in starting_rules:
			#d. Impose the restrictions we get from the other conjugate pairs
			for otheraut1, otheraut2 in system:
				extra_rule = rule.deduce(otheraut1, otheraut2)
				self.add(extra_rule)
	
	def iterate(self, system, message=None):
		if not self.atomic:
			self.atomise()
		old_count = len(self)
		self.deduce_new(system)
		if message is not None:
			message = "\n" + message
			print(message, self)
		self.atomise()
		if message is not None:
			print(message, self)
		if len(self) == old_count:
			warn("This iteration didn't produce new atoms")

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
		self.check_conjugate_pairs()
		self.list1_original = list1.copy()
		self.conjugator = Automorphism.identity((2,1))
		self.length = len(list1)
		
		# TODO: check for pairs (d_i, e_i) and (d_j, e_j) with either entry the same
		# if they are the same pair then remove one copy;
		# if they are different pairs the system is not SC
	
	@classmethod
	def from_dir(cls, dir='.', dname="d_{}.aut", ename="e_{}.aut", initial=0):
		#1. Determine the number of automorphism pairs in given directory
		length = 0
		while True:
			d = join(dir, dname.format(length + initial))
			e = join(dir, ename.format(length + initial))
			if not(exists(d) and exists(e)):
				break
			length += 1
		
		if length == 0:
			raise FileNotFoundError("Could not find at least one of {} or {} in '{}'".format(
				dname.format(initial), ename.format(initial), dir))
		list1 = [
			Automorphism.from_file( join(dir, dname.format(i + initial)) )
		for i in range(length)]
		list2 = [
			Automorphism.from_file( join(dir, ename.format(i + initial)) )
		for i in range(length)]
		
		return cls(list1, list2)
	
	@classmethod
	def random_solvable(cls, length=None):
		if length is None:
			length = randint(2, 5)
		list1 = [random_automorphism(signature=(2,1)) for _ in range(length)]
		conjugator = random_automorphism(signature=(2,1))
		list2 = [x^conjugator for x in list1]
		return cls(list1, list2)
	
	@classmethod
	def random(cls, length=None):
		if length is None:
			length = randint(2, 5)
		list1 = [random_automorphism(signature=(2,1))   for _ in range(length)]
		list2 = [x^random_automorphism(signature=(2,1)) for x in list1]
		return cls(list1, list2)
	
	def save_to_dir(self, dir='.'):
		for i, (d, e) in enumerate(self):
			d.save_to_file( join(dir, 'd_{}.aut'.format(i)) )
			e.save_to_file( join(dir, 'e_{}.aut'.format(i)) )
	
	def __len__(self):
		return self.length
	
	def __str__(self):
		return "<System of {} conjugacy equations>".format(self.length)
	
	def __getitem__(self, key):
		return self.list1[key], self.list2[key]
	
	def __iter__(self):
		yield from zip(self.list1, self.list2)
	
	def solved(self):
		return all(aut1 == aut2 for aut1, aut2 in self)
	
	def verify(self):
		return all(aut1^self.conjugator == aut2 for aut1, aut2 in zip(self.list1_original, self.list2))
	
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
	
	def check_conjugate_pairs(self):
		for index, (aut1, aut2) in enumerate(self):
			if not aut1.is_conjugate_to(aut2):
				raise ValueError("Pair #{} is not conjguate".format(index))
	
	def periodic_constraints(self):
		#Assert that period M orbits for aut1 have to go to period M orbits for aut 2
		for index, (aut1, aut2) in enumerate(self):
			for period in aut1.periodic_orbits:
				#a. Make a rule whose source is (aut1's period m elements)
				rule = Rule.empty()
				for orbit in aut1.periodic_orbits[period]:
					rule.source.extend(orbit)
				#b. The target for the rule is (aut2's periodic elements)
				for orbit in aut2.periodic_orbits[period]:
					rule.target.extend(orbit)
				rule.simplify()
				yield rule
	
	def align_periodic_orbits(self, iterations=1, verbose=False):
		rules = RuleSet(self.periodic_constraints())
		if verbose:
			print("Initial rules:", rules)
		if len(rules) == 0:
			warn("No rules deduced from periodic constraints")
			return rules
		for i in range(1, iterations+1):
			message = "Iteration {}".format(i) if verbose else None
			rules.iterate(self, message)
		return rules
	
	def coarsely_conjugate(self, atoms):
		if not atoms.atomic:
			raise ValueError("Given rules have not been atomised")
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
	
	def reset(self):
		self.list1 = self.list1_original
		self.conjugator = Automorphism.identity((2,1))

def _aut_in_V(x):
	return isinstance(x, Automorphism) and x.signature == (2, 1)

def bits(num, num_bits):
	for _ in range(num_bits):
		yield num % 2
		num >>= 1
