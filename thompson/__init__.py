"""Python scripts to work with elements of Thompson's groups F, T and V.

Importing ``*`` from this module will add the :class:`~thompson.word.Word`, :class:`~thompson.generators.Generators`, and :class:`~thompson.automorphism.Automorphism` classes to the current namespace. The :mod:`~thompson.examples` module is also imported.

.. moduleauthor:: David Robertson <david.m.robertson1@gmail.com>
"""

from .word import Word
from .generators import Generators
from .automorphism import Automorphism
from . import examples

#GLOBAL TODO
#Read through documentation and ensure that useful error messages are being raised, including for assertions.
