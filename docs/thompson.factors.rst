Free Factors and conjugacy
==========================

To :meth:`test for conjugacy <thompson.automorphism.Automorphism.test_conjugate_to>` we need to extract periodic and infinite components of an automorphism and test those components for conjugacy. These next few classes keep track of
- how the components embed into the original automorphism, and
- any extra information that we need to test the components for conjugacy. 

.. automodule:: thompson.factors

Factor class
------------

.. autoclass:: AutomorphismFactor
    :members:
    :undoc-members:

Periodic Factors
----------------

.. autoclass:: PeriodicFactor
    :members:
    :undoc-members:

Infinite Factors
----------------

.. autoclass:: InfiniteFactor
    :members:
    :undoc-members:

.. autofunction:: restricted_injections