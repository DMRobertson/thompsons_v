"""Python scripts to work with elements of Thompson's groups F, T and V.

Importing ``*`` from this module will add the :class:`~thompson.word.Word`, :class:`~thompson.generators.Generators`, and :class:`~thompson.automorphism.Automorphism` classes to the current namespace. The :mod:`~thompson.examples` module is also imported.

.. moduleauthor:: David Robertson <david.m.robertson1@gmail.com>
"""

from .word import Word
from .generators import Generators
from .automorphism import Automorphism
from .examples import *

"""
Global TODO:
- Read through documentation and ensure that useful error messages are being raised, including for assertions.
- Remove any unneccessary assertions?
- Remove from * imports where it's just me being lazy
- Change <A> to A* where neccessary. Need to check with AJD
- Having to fix quite a few 'forgot to copy' bug with mutable stuff (generating sets). Maybe copy everything at the start?

- Script to generate examples

"""
