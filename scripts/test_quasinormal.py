from test import setup_script
setup_script(__file__)

from thompson.word import *
from thompson.generators import Generators
from thompson.automorphism import Automorphism
from thompson.examples import *

u = Word('x', 2, 1)
example_4_25[u]

#Now let's see what we've computed
from pprint import pprint
d = example_4_25._dict
keys = list(sorted(d))
values = [d[key] for key in keys]

for key, value in zip(keys, values):
	print("{: <16} -> {:}".format(str(key), str(value)))