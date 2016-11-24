"""This is the top level of the package. Importing ``*`` from this module will add the following classes and functions to the current namespace.
"""

from .word            import Word
from .generators      import Generators
from .automorphism    import Automorphism
from .examples        import load_example, available_examples, show_examples, standard_generator
from .examples.random import random_automorphism
from .drawing         import plot, forest, flow, in_ipynb
from .cantorsubset    import CantorSubset
from .scsystem        import SCSystem

__all__ = [
	"Word",
	"Generators",
	"Automorphism",
	"load_example", "available_examples", "show_examples", "standard_generator",
	"random_automorphism",
	"plot", "forest", "flow",
	"SCSystem", "CantorSubset"
]

if in_ipynb():
	from IPython.display import display
	__all__.append("display")