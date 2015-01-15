MixedAuts
=============

.. automodule :: thompson.mixed

The MixedAut class
----------------------

.. autoclass:: thompson.mixed.MixedAut
    :members:
    :undoc-members:
    :exclude-members: inverse

.. todo::
	
	- Decide if the automorphism is in (the equivalent of) F, T, or V.
	- Add the named generators A, B, C, X_n of Thompson's V to the examples module. Analogues for general G_n,r?

Next steps
----------

To complete the details of the :meth:`conjugacy test <MixedAut.test_conjugate_to>`, we have to be able to test if two pure periodic automorphisms are conjugate and if two pure infinite automorphisms are conjugate. Any given automorphism can be broken down into a pure periodic and pure infinite part. These parts are called :meth:`free factors <MixedAut.free_factor>`.