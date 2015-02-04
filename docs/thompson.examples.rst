Examples
========

This module provides a number of explicit examples for use in doctests. Additionally, functions to generate random automorphisms are provided.

Explict named examples
----------------------

A list of named examples is loaded from the ``.aut`` files in `the examples folder <https://github.com/DMRobertson/thompsons_v/tree/master/thompson/examples>`_. This includes all the examples given in AJD and NB's paper.

.. todo:: Make Sphinx show these examples here? A better system for importing/loading examples?

.. automodule:: thompson.examples
    :members:
    :undoc-members:
    :show-inheritance:
    :exclude-members: random_signature, random_simple_word, random_basis, random_automorphism, random_periodic_automorphism, random_infinite_automorphism, random_conjugate_pair,  random_conjugate_periodic_pair, random_conjugate_infinite_pair, random_powers, random_power_conjugate_pair

Randomly generated examples
---------------------------

Words
^^^^^

.. autofunction:: random_signature
.. autofunction:: random_simple_word
.. autofunction:: random_basis

Automorphisms
^^^^^^^^^^^^^

.. autofunction:: random_automorphism
.. autofunction:: random_periodic_automorphism
.. autofunction:: random_infinite_automorphism
.. autofunction:: random_conjugate_pair
.. autofunction:: random_conjugate_periodic_pair
.. autofunction:: random_conjugate_infinite_pair
.. autofunction:: random_power_conjugate_pair
