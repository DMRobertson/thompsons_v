Automorphisms
=============

As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.

We consider three different classes of automorphism; this module implements functionality common to each class.

.. automodule :: thompson.automorphism

The Automorphisms class
----------------------

.. autoclass:: Automorphism
    :members:
    :undoc-members:
    :exclude-members: inverse

.. todo::
	
	- Decide if the automorphism is in (the equivalent of) F, T, or V.
	- Add the named generators A, B, C, X_n of Thompson's V to the examples module. Analogues for general G_n,r?

Next steps
----------

To complete the details of the :meth:`conjugacy test <MixedAut.test_conjugate_to>`, we have to be able to test if two pure periodic automorphisms are conjugate and if two pure infinite automorphisms are conjugate. Any given automorphism can be broken down into a pure periodic and pure infinite part. These parts are called :meth:`free factors <MixedAut.free_factor>`.