from fractions import Fraction

from .examples.random import random_automorphism
from .utilities       import is_dyadic

def test_conjugate_in_F(self, other):
	assert self.preserves_order()
	assert other.preserves_order()
	#1. Can the fixed point sets be match up?
	s_boundary = self.fixed_point_boundary()
	o_boundary = other.fixed_point_boundary()
	
	approx_conjugator = F_matching(s_boundary, o_boundary)
	if approx_conjugator is None:
		return None
	self = self ^ approx_conjugator
	if self.fixed_points() != other.fixed_points():
		return None
	
	#2. Apply the staircase algorithm

def F_matching(source, target):
	assert len(source) == len(target)
	assert source[0] == target[0] == 0
	assert source[-1] == target[-1] == 1
	
	last_source = Fraction(0, 1)
	last_target = Fraction(0, 1)
	domain = []
	range  = []
	
	for i in range(1, len(source)):
		next_source = source[i]
		next_target = target[i]
		if is_dyadic(next_source) != is_dyadic(next_target):
			return None
		if is_dyadic(next_source):
			subdomain = interval_to_words(last_source, next_source))
			subrange  = interval_to_words(last_target, next_target))
			subdomain.expand_to_size(len(subrange))
			subrange.expand_to_size(len(subdomain))
			domain.extend(subdomain)
			range.extend(subrange)
		else:
			pass
			#use Prop 4.16 to build an element as required
		last_source = next_source
		last_target = next_target
	
	domain = Generators((2, 1), domain)
	range  = Generators((2, 1), range)
		
	return Automorphism(domain, range)

def interval_to_words(left, right):
	assert is_dyadic(left)
	assert is_dyadic(right)
	assert left < right
	leftmost = dyadic_as_side(left, LHS=True)
	rightmost = dyadic_as_side(right, LHS=False)
	yield from words_between(leftmost, rightmost)
	
def dyadic_as_side(target, LHS=True):
	"""Express the dyadic value *target* as the LHS or RHS of some dyadic interval [a/2^m, (a+1)/2^m] in the form of a Word."""
	assert is_dyadic(target)
	word  = Word("x", (2, 1)) # word with endpoints [left, right]
	left  = Fraction(0, 1)
	mid   = Fraction(1, 2)
	right = Fraction(1, 1)
	current = None
	while current != target:
		if target < mid or (target == mid and not LHS):
			word = word.expand(0)
			end = mid
		else:
			word = word.expand(1)
			start = mid
		mid = (left + right) / 2;
		current = left if LHS else right
	return word

def words_between(leftmost, rightmost):
	assert leftmost <= rightmost
	word = leftmost
	...
	#yield ...

def main():
	g = random_automorphism((2, 1), group='F')
	h = random_automorphism((2, 1), group='F')
	k = g ^ h
	test_conjugate_in_F(g, k)