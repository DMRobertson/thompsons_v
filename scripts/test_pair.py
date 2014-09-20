import sys, os
sys.path.insert(0, os.path.abspath('.'))

from thompson.tree_pair import TreePair
from thompson.trees import DrawableTree

#Testing a composition of tree pairs
f = TreePair("100", "100", "2 1")
f.range.name = "A"
f.range.left.name = "B"
f.range.right.name = "C"

f.domain.name = "A'"
f.domain.left.name = "C'"
f.domain.right.name = "B'"

g = TreePair("10100", "10100", "3 1 2")
g.domain.name = "F"
g.domain.left.name = "G"
g.domain.right.name = "H"
g.domain.right.left.name = "I"
g.domain.right.right.name = "J"

g.range.name = "F'"
g.range.left.name = "J'"
g.range.right.name = "H'"
g.range.right.left.name = "G'"
g.range.right.right.name = "I'"

print(f*g == TreePair("11000", "10100", "3 2 1"))