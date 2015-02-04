Free Factors and conjugacy
==========================

To :meth:`test for conjugacy <thompson.mixed.MixedAut.test_conjugate_to>` we need to extract periodic and infinite components of an automorphism and test those components for conjugacy. These next few classes keep track of

- how the components embed into the original automorphism, and
- any extra information that we need to test the components for conjugacy. 

Periodic Factors
----------------

.. automodule:: thompson.periodic
    :members:
    :undoc-members:

Infinite Factors
----------------

.. automodule:: thompson.infinite
    :members:
    :undoc-members:

Next steps
----------

With these classes, the conjugacy and power conjugacy tests are implemented. The other important part of this package is the :mod:`~thompson.examples` module.