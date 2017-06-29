def pairwise(iterable):
	#Nicked from itertools
	"s -> (s0,s1), (s1,s2), (s2, s3), ..."
	a, b = tee(iterable)
	next(b, None)
	return zip(a, b)

def increasing_sequence(iterable):
	return all(
		prev < next
		for prev, next in pairwise(iterable)
	)

def all_satisfy(iterable, func):
	return all( func(x) for x in iterable )

def grad(start, end):
	return (end[1] - start[1]) / (end[0] - start[0])
	
def int_power_of_two(n):
	#Assume n >= 1
	return n & (n - 1) == 0