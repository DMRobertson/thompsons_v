Automorphisms
=============

As usual, a bijective homomorphism is called an *isomorphism*, and an isomorphism from an algebra to itself is called an *automorphism*. A collection of automorphisms of an object :math:`O` forms a group :math:`\mathrm{Aut}(O)` under composition. The group of automorphisms :math:`\mathrm{Aut}(V_{n,r})` is known as :math:`G_{n,r}`.

We consider three different classes of automorphism; this module implements functionality common to each class.

.. automodule :: thompson.automorphism

The Automorphisms class
-----------------------

.. autoclass:: Automorphism
    :members:
    :undoc-members:
    :exclude-members: inverse

Membership placeholders
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: thompson.membership

Next steps
----------

.. currentmodule:: thompson.automorphism

Because an Automorphism can be :meth:`repeatedly applied <Automorphism.__pow__>`, we may consider the orbit :math:`\{w\psi^n \mid n \in \mathbb{Z}\}` of any word :math:`w`. (Note that this is the orbit of :math:`w` under the cyclic subgroup :math:`\langle \psi \rangle`, rather than all of :math:`G_{n,r}`.)

The next module gives us tools to classify and work with these orbits. (In truth, the :class:`Automorphism` class uses these tools, so this documentation is really in the wrong order.)