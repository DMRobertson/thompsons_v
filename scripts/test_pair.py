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

x = TreePair("110100100", "110010100", "4 5 1 2 3")
print(x.in_F(), x.in_T())
x.render(filename='original.svg')

x.reduce()
x.render(filename='reduced.svg')

y = TreePair("100", "100", "2 1")
y.render(filename='expected.svg')

TreePair("11000", "10100", "1 2 3").render(filename='test.svg')
TreePair("11000", "10100", "1 2 3").render_bijection(filename='bijection.svg')