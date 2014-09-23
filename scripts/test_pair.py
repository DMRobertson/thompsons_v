import sys, os
sys.path.insert(0, os.path.abspath('.'))

from thompson.tree_pair import TreePair

a = TreePair("11010100100", "11010100100", "1 6 5 2 3 4")
x = a
x.render(filename='a_power_1.svg')
for i in range(5):
	x = x * a
	x.render(filename='a_power_{}.svg'.format(i+2))

#TODO: __pow__: method of squares?
y = a * a * a * a * a * a
print('a**6 == 1?', x == TreePair("0", "0"), y == TreePair("0", "0"))
