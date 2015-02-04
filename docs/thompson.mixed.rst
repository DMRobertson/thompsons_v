MixedAuts
=============

.. automodule :: thompson.mixed

A general automorphism :math:`\phi \in G_{n,r}` may exhibit periodic behaviour, infinite behaviour, or a mix of both. We say that an automorphism is

- *purely periodic* if it has finite order;
- *purely infinite* if it has no periodic :mod:`~thompson.orbits`; and
- *mixed* otherwise.

We can represent a mixed automorphism :math:`\psi` as a `free product <http://en.wikipedia.org/wiki/Free_product>`_ :math:`\psi = \psi_P * \psi_I` of a purely periodic and purely infinite automorphism. The automorphisms :math:`\psi_P` and :math:`\psi_I` are called :meth:`~MixedAut.free_factors`. We form this decomposition because the conjugacy problem is easier to solve when the automorphisms are both pure infinite or both pure periodic.

This class is responsible for:

- Generating the free factors from a mixed automorphism;
- Combining a periodic and infinite factor into a mixed automorphism; and
- Delegating the conjugacy and power conjugacy tests to these factors.

The MixedAut class
------------------

.. autoclass:: thompson.mixed.MixedAut
    :members:
    :undoc-members:
    :exclude-members: inverse

Next steps
----------

To complete the details of the :meth:`conjugacy test <MixedAut.test_conjugate_to>`, we have to be able to test if two pure periodic automorphisms are conjugate and if two pure infinite automorphisms are conjugate. Any given automorphism can be broken down into a pure periodic and pure infinite part. These parts are called :meth:`free factors <MixedAut.free_factor>`.