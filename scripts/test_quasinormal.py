from test import setup_script
setup_script(__file__)

from thompson.word import Word, from_string, standardise, format
from thompson.generators import Generators
from thompson.automorphism import Automorphism
# from thompson.examples import *

#Extraction only
print('word is: ', Word("x a1 a2 a1 x a2 a1 L a1 a2", 2, 1))
print('expected:', 'x1 a1 a2 a1 a2')
print()
#Contraction only
print('word is: ', Word("x2 a2 a2 a1 x2 a2 a2 a2 L", 2, 2))
print('expected:','x2 a2 a2')
print()
#Contraction only, arity 3 
print('word is: ', Word("x1 a1 a1 x1 a1 a2 x1 a1 a3 L", 3, 2))
print('expected:', 'x1 a1')
print()
"""
print(Word("x a1 a2 a1 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L", 2, 1))
print('x1 a1 a2 a1')
print('='*5)
print()

print(Word("x a1 a2 a1 x a1 a2 a1 x a2 a1 L a1 a2 a1 x a1 a2 a1 a2 a2 L L", 2, 1))
print('x1 a1 a2 a1 x1 a1 a2 a1 a2 L')
print('='*5)
print()

print(Word("x a2 x a1 L x a1 L a1 a2 a1 a2", 2, 1))
print('x1 a1 a1 a2')
print('='*5)
print()
"""

# w = Word('x1 a1 a2', 2, 1)
# print(example_4_25[w])