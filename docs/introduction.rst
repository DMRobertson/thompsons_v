Getting Started
===============

Installation
------------

1. Python interpreter.
	This program is written in `Python <python.org>`_ and needs a Python interpreter to run; to be specific, **Python 3.3 or above** is required. If you don't have this already installed on your system, your best bet is to download the latest version (currently 3.5) from the `Python website <python.org>`_. For Linux users, Python may be available from your distribution's repositories; take care to install version 3.x rather than 2.x!

2. Source code.
	If you have `git <git-scm.com>`_ available on your system, the easiest way to obtain the program is to simply clone the repository:
		
		``git clone https://github.com/DMRobertson/thompsons_v.git``
		
	If this is not possible, a ZIPped version can be downloaded from the `project page <https://github.com/DMRobertson/thompsons_v>`_ on GitHub.

3. Package requirements.
	The de-facto package management tool for Python is the ``pip`` program, which is installed as part of Python 3.4 and higher. Python 3.3 users should consult the `pip documentation <https://pypi.python.org/pypi/pip/>`_ for instructions on how to install ``pip`` directly. Some Linux distributions may also provide pip as a package (make sure it's for Python 3---it may be listed as ``pip3``).
	
	Once pip is installed, navigate to the project directory and run
	
		``pip install -r requirements.txt``
	
	(possibly you may need to use ``pip3``). This should do all the hard work for you.
	
	Direct installation:
		If you would prefer not to install ``pip``, you may install the prerequisites directly. There is only one major Python package required to run the program: `NetworkX <https://networkx.github.io/>`_, a library for working with graphs in Python. This may be available in Linux repositories, e.g. `Ubuntu <http://packages.ubuntu.com/search?keywords=python3-networkx&searchon=names>`_.
		
		This documentation is built using `Sphinx <http://sphinx-doc.org/>`_. If you wish to build it yourself, or run the test suite, you must install Sphinx also. Again, consider checking Linux repositories e.g. `in Ubuntu <http://packages.ubuntu.com/search?suite=default&section=all&arch=any&keywords=python3-sphinx&searchon=names>`_.

First steps
-----------

Navigate to the root folder (``thompsons_v``) in a terminal and start up Python using the ``python`` command. Once the python prompt (``>>>``) is available, try to ``import thompson``::

	$ python
	Python 3.3.3 (v3.3.3:c3896275c0f6, Nov 18 2013, 21:19:30) [MSC v.1600 64 bit (AMD64)] on win32
	Type "help", "copyright", "credits" or "license" for more information.
	>>> import thompson
	>>>

If this fails, check that you're using Python 3.3 or above. If ``python`` started up Python 2.x, try using the command ``python3`` instead.

We can now use the python `REPL <http://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop>`_ as a sort of calculator for these groups.

.. doctest::
	:options: +NORMALIZE_WHITESPACE
	
	>>> from thompson import *
	>>> sig = (2, 1) #signature of the algebra we work in
	>>> domain = Generators.standard_basis(sig).expand(0).expand(0)
	>>> range  = Generators.standard_basis(sig).expand(0).expand(1)
	>>> print(domain)
	[x1 a1 a1, x1 a1 a2, x1 a2]
	>>> print(range)
	[x1 a1, x1 a2 a1, x1 a2 a2]
	>>> phi = Automorphism(domain, range)
	>>> print(phi)
	InfiniteAut: V(2, 1) -> V(2, 1) specified by 3 generators (after expansion and reduction).
	x1 a1 a1 -> x1 a1   
	x1 a1 a2 -> x1 a2 a1
	x1 a2    -> x1 a2 a2

.. doctest::
	:options: +NORMALIZE_WHITESPACE
	
	>>> phi.image('x a2 a1 a2 a2 a1')
	Word('x1 a2 a2 a1 a2 a2 a1', (2, 1))
	>>> phi.order
	inf
	>>> phi.print_characteristics()
	(-1, a1)
	(1, a2)
	>>> phi.dump_QNB()
	x1 a1 Left semi-infinite component with characteristic (-1, a1)
	x1 a2 Right semi-infinite component with characteristic (1, a2)

.. doctest::
	:options: +NORMALIZE_WHITESPACE
	
	>>> print(phi ** 2)
	InfiniteAut: V(2, 1) -> V(2, 1) specified by 4 generators (after expansion and reduction).
	x1 a1 a1 a1 -> x1 a1      
	x1 a1 a1 a2 -> x1 a2 a1   
	x1 a1 a2    -> x1 a2 a2 a1
	x1 a2       -> x1 a2 a2 a2
	>>> print(~phi) #inverse
	InfiniteAut: V(2, 1) -> V(2, 1) specified by 3 generators (after expansion and reduction).
	x1 a1    -> x1 a1 a1
	x1 a2 a1 -> x1 a1 a2
	x1 a2 a2 -> x1 a2   
	>>> print(~phi * phi) # == identity
	PeriodicAut: V(2, 1) -> V(2, 1) specified by 1 generators (after expansion and reduction).
	x1 -> x1

Locating examples
-----------------

A :mod:`number of examples <thompson.examples>` are included in this package, some of which are used in [BDR]_ . To access them, use the :func:`~thompson.examples.load_example` function:

.. doctest::
	:options: +NORMALIZE_WHITESPACE
	
	>>> from thompson import *
	>>> phi = load_example('nathan_pond_example')
	>>> print(phi)
	InfiniteAut: V(2, 1) -> V(2, 1) specified by 7 generators (after expansion and reduction).
	x1 a1 a1 a1 a1 -> x1 a1 a1      
	x1 a1 a1 a1 a2 -> x1 a1 a2 a1 a1
	x1 a1 a1 a2 a1 -> x1 a2 a2      
	x1 a1 a1 a2 a2 -> x1 a2 a1      
	x1 a1 a2       -> x1 a1 a2 a1 a2
	x1 a2 a1       -> x1 a1 a2 a2 a2
	x1 a2 a2       -> x1 a1 a2 a2 a1

.. note:: Previously, one would access this example by using ``from thompson.examples import nathan_pond_example``. I changed this, because it meant that every example was loaded (and the quasinormal bases etc. computed) whenever ``thompson.examples`` was imported.

To see the list of available examples, consult the :mod:`documentation for the examples module <thompson.examples>`. You might also like to consult `the demo notebook <https://github.com/DMRobertson/thompsons_v/blob/master/notebooks/Thompson%20example.ipynb>`_ on GitHub.


Running the test suite
----------------------

Make sure Sphinx is installed (see the installation section). Then navigate to the ``docs`` folder and run ``make doctest``. This should look through the source code for short tests and alert you if any of them fail. You can also run ``make html`` to build this documentation. Reset with ``make clean``.
