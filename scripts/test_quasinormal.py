from test import setup_script
setup_script(__file__)

from thompson.word import Word
from thompson.examples import *

example_4_5.dump_mapping()
orbit_types(example_4_5, example_4_5.domain)
