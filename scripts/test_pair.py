#Naughty hack to make the thompson package importable. Adds the root directory to sys.path
import sys
from pathlib import Path

root = Path(__file__).parents[1]
sys.path.insert(0, str(root))

from thompson.tree_pair import TreePair

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
