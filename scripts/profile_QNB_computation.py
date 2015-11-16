from scripts import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson import Automorphism
from thompson.examples import *

filename = 'QNB_size_206'
filename = 'QNB_size_344'
aut = Automorphism.from_file(filename + '.aut')

profile = 'cProfile'

print(aut)
#TODO. This is broken since I made the QNB automatically computed in __init__.
if profile == 'cProfile':
    import cProfile
    print('profiling with cProfile')
    cProfile.run('X = aut.quasinormal_basis')#, filename=filename+'.profile')
    
elif profile == 'callgraph':
    from pycallgraph import PyCallGraph
    from pycallgraph.output import GraphvizOutput
    
    print('profiling with PyCallGraph')
    with PyCallGraph(output=GraphvizOutput()):
        X = aut.quasinormal_basis
else:
    X = aut.quasinormal_basis

print(X)
for x in X:
    ctype, _, _ = aut.orbit_type(x, X)
    print(x, ctype)