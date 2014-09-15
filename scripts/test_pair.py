import sys, os
sys.path.insert(0, os.path.abspath('.'))

from thompson.tree_pair import TreePair
from thompson.trees import DrawableTree

#Testing a composition of tree pairs
#The same number of leaves, but not compatible
# f = TreePair("1101000", "1100100")
# g = TreePair("1110000", "1011000")
#want the tree pair f*g, corresponding to the bijection g after f of [0, 1]
# h = f * g
# print(h)

#1. result = new tree pair: domain = f.domain, range = g.range
#2. adjust so that they are compatible.
# for each leaf of f.range:
	#does a node with the same path (corresponding to the same interval) exist in g.domain?
	#if so, and the node is a leaf: all is well, continue
	#if so, but the node is not a leaf: duplicate the structure below the node in g.domain in f.range
		#and update permutation
	#if not: add children to g.domain 

"""
f = TreePair("1110010010100", "1101010011000")
f.render(filename='original.svg', name='f')
f.reduce()
f.render(filename='reduced.svg', name = 'f reduced')

TreePair("11100010100", "11010011000").render(filename='expected.svg', name='expected')
"""

# f = TreePair("111100000", "111100000", "1 2 3 5 4")
# f = TreePair("11000", "11000", "1 3 2")
"""
f = TreePair("11101000101100100", "11011000101010100", "3 7 6 1 2 4 5 8 9")
f.render(filename='original.svg', name='f')
f.reduce()
f.render(filename='reduced.svg', name = 'f reduced')
TreePair("111010001011000", "110110001010100", "3 7 6 1 2 4 5 8").render(filename='expected.svg', name='f expected')"""

f = TreePair("11101000101100100", "11011000101010100", "3 7 6 1 2 4 5 8 9")
f.render(filename='original.svg', name='g')
f.reduce()
f.render(filename='reduced.svg', name = 'g reduced')