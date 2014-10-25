Automorphisms and conjugacy
===========================

.. automodule :: thompson.automorphism

The Isomorphism class
---------------------

.. autoclass:: thompson.automorphism.Isomorphism
    :members:
    :undoc-members:

The Automorphism class
----------------------

.. autoclass:: thompson.automorphism.Automorphism
    :members:
    :undoc-members:

Periodic and infinite components
--------------------------------

To :meth:`test for conjugacy <Automorphism.test_conjugate_to>` we need to extract periodic and infinite components of an automorphism and test those components for conjugacy. These next few classes keep track of
- how the components embed into the original automorphism, and
- any extra information that we need to test the components for conjugacy. 

.. autoclass:: AutomorphismFactor
    :members:
    :undoc-members:

.. autoclass:: thompson.automorphism.PeriodicFactor
    :members:
    :undoc-members:

.. autoclass:: thompson.automorphism.InfiniteFactor
    :members:
    :undoc-members:

Next steps
----------

.. todo:: Introduce the orbits module here.