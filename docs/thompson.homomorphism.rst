Homomorphisms
=============

.. automodule :: thompson.homomorphism
    :members:
    :undoc-members:
    :exclude-members: Homomorphism

The Homomorphism Class
----------------------

.. autoclass :: thompson.homomorphism.Homomorphism
    :members:
    :undoc-members:


.. rubric:: **Footnotes**

.. [#footnote_why_optional_reduce] Sometimes we have to expand :class:`free factors <thompson.factors.AutomorphismFactor>` so that the orbit sizes match up.

.. [#footnote_why_optional_image_args] Exposing more arguments means I can write this function in a more general manner. Doing so makes it easy to compute :meth:`inverse images <thompson.automorphism.Automorphism.image>` under an :class:`~thompson.automorphism.Automorphism`.

Next steps
----------

It is not possible to repeatedly apply a homomorphism :math:`\psi` to a word unless :math:`\psi` is actually an automorphism.
Then we can consider the orbit :math:`\{w\psi^n \mid n \in \mathbb{Z}\}` of any word :math:`w`. The next module gives us tools to classify and work with these orbits.
