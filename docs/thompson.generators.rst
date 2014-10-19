Generating sets and bases 
=========================

.. automodule:: thompson.generators
    :members:
    :undoc-members:
    :show-inheritance:

Next steps
----------

Suppose we have a map :math:`f \colon V_{n,r} \to V_{n',r'}`. If this is just a set map then we would have to specify the image of every word in :math:`V_{n,r}`. However, if :math:`f` preserves structure then it is sufficient to specify the image of a generating set---or even better, a basis---for :math:`V_{n,r}`.

The next module describes automorphisms as bijections between bases. This is the bulk of the package, and implements the algorithms described in the paper.