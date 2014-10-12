from test import setup_script
dwg, canvas = setup_script(__file__)

from thompson.tree_pair import *
from thompson.drawing import *

id = TreePair("0", "0")

a = TreePair("11010100100", "11010100100", "1 6 5 2 3 4")
x = a
x.render(filename='a_power_1.svg')
for i in range(5):
	x = x * a
	x.render(filename='a_power_{}.svg'.format(i+2))
	print('a to the power {} is the identity? {}'.format(i+2, x == id))

#TODO: __pow__: method of squares?
#Does one big chain of multiplications work?
y = a * a * a * a * a * a
print('a**6 == 1?', y == id)