"""This is the top level of the package. Importing ``*`` from this module will add the following classes and functions to the current namespace.
"""

from .word            import Word
from .generators      import Generators
from .automorphism    import Automorphism
from .examples        import load_example, available_examples, show_examples
from .examples.random import random_automorphism

__all__ = ["Word", "Generators", "Automorphism", "load_example", "available_examples", "show_examples", "random_automorphism"]