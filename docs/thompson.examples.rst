Examples
========

This module provides a number of explicit examples for use in doctests. Additionally, functions to generate random automorphisms are provided.

Explict named examples
----------------------

A list of named examples is loaded from the ``.aut`` files in `the examples folder <https://github.com/DMRobertson/thompsons_v/tree/master/thompson/examples>`_. This includes all the examples given in AJD and NB's paper.

.. todo:: Make Sphinx show these examples here?

.. automodule:: thompson.examples
    :members:
    :undoc-members:
    :show-inheritance:
    :exclude-members: random_signature, random_simple_word, random_basis, random_automorphism, random_conjugate_pair, random_conjugate_factors, random_conjugate_periodic_factors, random_conjugate_infinite_factors

Randomly generated examples
---------------------------

.. autofunction:: random_signature
.. autofunction:: random_simple_word
.. autofunction:: random_basis
.. autofunction:: random_automorphism
.. autofunction:: random_conjugate_pair
.. autofunction:: random_conjugate_factors
.. autofunction:: random_conjugate_periodic_factors
.. autofunction:: random_conjugate_infinite_factors