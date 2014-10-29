.. currentmodule:: thompson.word

Words and standard forms
========================

By a *word*, we mean a string written using the symbols

.. math:: X \cup \Omega = \{x_1, \dotsc, x_r\} \cup \{\alpha_1, \dotsc, \alpha_n, \lambda\}.

The :math:`x_i` are constants, the :math:`\alpha_i` are unary operators and :math:`\lambda` is an *n*-ary operation. In this package, we refer to the parameters :math:`n` and :math:`r` as the *arity* and *alphabet size* respectively.

If we collect together all such words, we obtain an *algebra* (in the sense of `universal algebra <http://en.wikipedia.org/wiki/Universal_algebra>`_) of words. We consider three different algebras; using the notation of [Cohn]_:

1. :math:`W(\Omega; X)`, the set of all finite strings written using letters in :math:`X \cup \Omega`. Cohn calls these :math:`\Omega`-rows.
2. :math:`W_\Omega(X)`, the subset of :math:`W(\Omega; X)` whose strings begin with an :math:`x_i`, and represent a :func:`valid <validate>` series of operations. Cohn calls these :math:`\Omega`-words.
3. :math:`V_{n, r}(X) = W_\Omega(X) / \mathfrak{q}`. This is equivalent to the set of words in :math:`W_\Omega(X)` which are in Higman's :py:func:`standard form <thompson.word.standardise>`. Roughly speaking, a word is in standard form if is valid (type 2) and it cannot be reduced to a shorter word using the following two rules. 
	
	a. :math:`w_1 \dots w_n \lambda \alpha_i = w_i`, where each :math:`w_i` is in standard form.
	b. :math:`w \alpha_1 \dots w \alpha_n \lambda = w`, where :math:`w` is in standard form.

Sometimes we need to refer to a string which consists only of :math:`\alpha`-s. Write :math:`A` for the set :math:`\{\alpha_1, \dotsc, \alpha_n\}`. We define :math:`\langle A \rangle` to be the set of finite strings over :math:`A`.

.. seealso:: Definition 2.1, Section 2.3, Remark 3.3 and Definition 3.7 of the paper.

Signatures
----------

To each word :math:`w \in V_{n,r}` we associate a :class:`Signature` object. This allows us to pass words around as arguments to functions and keep track of which algebra :math:`V_{n,r}` the word belongs to. Signatures are implemented as glorified :func:`namedtuples <py3:collections.namedtuple>`. In particular this means they are `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_.

.. autoclass:: thompson.word.Signature
    :members:
    :undoc-members:

Operations on words
-------------------

We represent a word as a :class:`tuple <py3:tuple>` of integers, where:

- :math:`x_i` is represented by :math:`i`,
- :math:`\alpha_i` is represented by :math:`-i`, and
- :math:`\lambda` is represented by :math:`0`.

We can write words of all types in this format, but we're only interested in the standard forms (type 3). To this end, the :func:`~validate` detects those which are of type 2, and :func:`~standardise` reduces type 2 words to type 3.

.. automodule:: thompson.word
    :members:
    :undoc-members:
    :exclude-members: Word, Signature

The Word class
--------------

It's important to know when a sequence of letters denotes a word in standard form. The :class:`Word` class addresses this problem by standardising its input upon creation. We can think of a Word object as a fixed list of letters with the guarantee that its are in standard form. Words are also given a :class:`Signature` at creation time, so that we know which algebra the word comes from.

We also need to know that once a Word has been created, it cannot be modified to something not in standard form. We achieve this simply by making it impossible to modify a Word. Words are implemented as (subclasses of) :class:`tuple <py3:tuple>`, which are `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_. One useful side effect of this is that Words can be used as dictionary keys. 

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

Next steps
----------

Now that we can work with words in :math:`V_{n,r}`, we need to be able to work with *collections* of words. The :mod:`~thompson.generators` module will let us treat a list of words as a generating set and examine the subalgebra that the collection generates.