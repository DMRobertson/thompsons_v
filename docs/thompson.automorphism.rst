Automorphisms
=============

.. automodule :: thompson.automorphism

The Automorphism class
----------------------

.. autoclass:: thompson.automorphism.Automorphism
    :members:
    :undoc-members:
    :exclude-members: inverse


Next steps
----------

To complete the details of the :meth:`conjugacy test <Automorphism.test_conjugate_to>`, we have to be able to test if two pure periodic automorphisms are conjugate and if two pure infinite automorphisms are conjugate. Any given automorphism can be broken down into a pure periodic and pure infinite part. These parts are called *free factors*.