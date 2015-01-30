.. currentmodule:: thompson.word

Words and standard forms
========================

By a *word*, we mean a `string <http://en.wikipedia.org/wiki/Formal_language#Words_over_an_alphabet>`_ written using the symbols

.. math:: X \cup \Omega = \{x_1, \dotsc, x_r\} \cup \{\alpha_1, \dotsc, \alpha_n, \lambda\}.

We treat the :math:`x_i` are constants, the :math:`\alpha_i` are unary operators and :math:`\lambda` as an *n*-ary operation. We refer to the parameters :math:`n` and :math:`r` as the *arity* and *alphabet size* respectively.

Notation
--------

If we collect together all such words, we obtain an *algebra* (in the sense of `universal algebra <http://en.wikipedia.org/wiki/Universal_algebra>`_) of words. We consider three different algebras; using the notation of [Cohn]_:

.. _type1:

1. :math:`W(\Omega; X)`, the set of all finite strings written using letters in :math:`X \cup \Omega`. In other words/jargon, this is the `free monoid <http://en.wikipedia.org/wiki/Free_monoid>`_ on :math:`X \cup \Omega`. Cohn calls these strings :math:`\Omega`-rows.

.. _type2:

2. :math:`W_\Omega(X)`, the subset of :math:`W(\Omega; X)` whose strings represent a :func:`valid <validate>` series of operations. *Valid* means that each the operations :math:`\alpha_i` and :math:`\lambda` always recieve the correct number of arguments. Cohn calls these strings :math:`\Omega`-words.

.. _type3:

3. :math:`V_{n, r}(X) = W_\Omega(X) / \mathfrak{q}`. This is equivalent to the set of words in :math:`W_\Omega(X)` which are in Higman's :py:func:`standard form <thompson.word.standardise>`. Sections 2 and 3 of the paper give the theory. Roughly speaking, a word is in standard form if is valid (type 2) and it cannot be reduced to a shorter word using the following two rules. 
	
	a. :math:`w_1 \dots w_n \lambda \alpha_i = w_i`, where each :math:`w_i` is in standard form.
	b. :math:`w \alpha_1 \dots w \alpha_n \lambda = w`, where :math:`w` is in standard form.

.. comment: We usually drop the :math:`X` in :math:`V_{n, r}(X)` when it is clear what the . This is because we don't really need to keep track of the set containing the elements :math:`x_1, \dotsc x_n` which are the constants of our algebra.

Sometimes we need to refer to a string which consists only of :math:`\alpha`-s. Write :math:`A` for the set :math:`\{\alpha_1, \dotsc, \alpha_n\}`. We define :math:`A^*` to be the set of finite strings over :math:`A`. (Formally this is the `free monoid <http://en.wikipedia.org/wiki/Free_monoid>`_ on :math:`A`.)

Finally, let :math:`S` be a set of words. We define :math:`X\langle A\rangle` to be the set of concatenations :math:`s \Gamma`, where :math:`s \in S` and :math:`\Gamma \in A^*`. It is sometimes helpful to think of :math:`S\langle A\rangle` as the set of words *below* :math:`S`.

.. seealso:: Definition 2.1, Section 2.3, Remark 3.3 and Definition 3.7 of the paper.

Signatures
----------

It may happen that we need to work in different algebras :math:`V_{n,r}, V_{m,s}` at the same time. To keep track of the algebras that different words belong to, we associate a :class:`Signature` to each word :math:`w \in V_{n,r}`. [#footnote_why_signatures]_ Signatures are implemented as a glorified :math:`(n, r)`. This means they are `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_.

.. autoclass:: thompson.word.Signature
    :members:
    :undoc-members:

Operations on words
-------------------

We represent a word as a :class:`tuple <py3:tuple>` of integers, where:

- :math:`x_i` is represented by :math:`i`,
- :math:`\alpha_i` is represented by :math:`-i`, and
- :math:`\lambda` is represented by :math:`0`.

We can write words of all types in this format, but we're only interested in the standard forms (:ref:`type 3 <type3>`). To this end, the :func:`~validate` detects those which are of :ref:`type 2 <type2>`, and :func:`~standardise` reduces type 2 words to type 3.

.. automodule:: thompson.word
    :members:
    :undoc-members:
    :exclude-members: Word, Signature

The Word class
--------------

It's important to know when a sequence of letters denotes a word in standard form. The :class:`Word` class addresses this problem by standardising its input upon creation. We can think of a Word object as a fixed list of letters with the guarantee that its are in standard form. Words are also given a :class:`Signature` at creation time, so that we know which algebra the word comes from.

We also need to know that once a Word has been created, it cannot be modified to something not in standard form. We achieve this simply by making it impossible to modify a Word. Words are implemented as (subclasses of) :class:`tuple <py3:tuple>`, which are `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_. One useful side effect of this is that Words can be used as dictionary keys. [#footnote_immutable_words]_

.. doctest::
	
	>>> w = Word('x1 a1', (2, 1))
	>>> x = {}; x[w] = 'stored value'
	>>> x[w]
	'stored value'
	>>> #Can also use the underlying tuple as a key
	>>> x[1, -1]
	'stored value'
	>>> #the reason this works
	>>> hash(w) == hash((1, -1))
	True

.. autoclass:: Word
    :members:
    :undoc-members:

.. rubric:: **Footnotes**

.. [#footnote_why_signatures] This made it easier to catch bugs before they happened when I was passing words around as arguments to functions. 
.. [#footnote_immutable_words] This made it so much easier to cache the images of :class:`homomorphisms <thompson.homomorphism.Homomorphism>`, because I could just use a vanilla Python dictionary.

Next steps
----------

Now that we can work with words in :math:`V_{n,r}`, we need to be able to work with *collections* of words. The :mod:`~thompson.generators` module will let us treat a list of words as a generating set and examine the subalgebra that the collection generates.