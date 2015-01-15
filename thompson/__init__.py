"""Python classes to work with the Higman-Thompson groups :math:`G_{n,r}`.

Importing ``*`` from this module will add the :class:`~thompson.word.Word`, :class:`~thompson.generators.Generators`, and :class:`~thompson.mixed.MixedAut` classes to the current namespace. The :mod:`~thompson.examples` module is **not** imported.
"""

from .word import Word
from .generators import Generators
from .mixed import MixedAut
