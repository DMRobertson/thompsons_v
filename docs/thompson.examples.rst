Examples
========

.. automodule:: thompson.examples

Explict named examples
----------------------

A list of named examples is loaded from the ``.aut`` files in `the examples folder <https://github.com/DMRobertson/thompsons_v/tree/master/thompson/examples>`_.
This includes all the examples given in the paper, as well as others used to test the package.
Use the following functions to load one of these examples.

.. autofunction:: available_examples
.. autofunction:: load_example
.. autofunction:: load_example_pair
.. autofunction:: load_all_examples
.. autofunction:: standard_generator

.. todo:: Add the named generators A, B, C of Thompson's T and V. Analogues for general G_n,r?

.. _examples-table:

List of named examples
----------------------

Note that some automorphisms have more than one name---they will appear once in this list for every alias.

.. todo:: make the documentation build process include links to forest pair diagrams for these elements.

.. include:: examples_table.txt

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
