A package for Python 3.3+ to work with elements of Thompson's groups F, T and V. Documentation is built automatically at [Read The Docs](http://thompsons-v.readthedocs.org/en/latest/).

![TreePair("11000", "10100", "1 2 3")](https://rawgit.com/DMRobertson/thompsons_v/master/docs/examples/tree_pair/TreePair_mul.svg "Multiplication in V.")


Requirements
------------

The package code itself only uses two packages:
>  - [svgwrite](https://pypi.python.org/pypi/svgwrite/)
>  - [pathlib](https://pypi.python.org/pypi/pathlib/) (included in Python 3.4's standard library).

If either package is missing, try running ``pip install <packagename>`` or ``pip3 install <packagename>`` to install locally.

The HTML documentation is generated using [Sphinx](https://pypi.python.org/pypi/Sphinx). To build, ensure that Sphinx is installed (``pip install sphinx``), then navigate to [the docs folder](/docs) and run ``make html``. If nothing happens, run ``make clean`` and try ``make html`` again.

More interestingly, I'm using Sphinx's [doctest extension](http://sphinx-doc.org/ext/doctest.html) to run some simple automated tests. To use it, ensure Sphinx is installed and run ``make doctest`` from the docs directory.


How do I start playing around?
------------------------------

Either clone the repository with

	git clone https://github.com/DMRobertson/thompsons_v.git

or download [the ZIPped version](https://github.com/DMRobertson/thompsons_v/archive/master.zip) and extract it. Navigate to the root directory---the one containing ``docs``, ``scripts`` and ``thompson``. Start the Python interpreter and run:

	>>> from thompson.tree_pair import TreePair
	>>> x = TreePair("11000", "10100", "1 2 3")
	>>> x.render(filename='example.svg')

This will save the following picture in the file ``example.svg``.

![TreePair("11000", "10100", "1 2 3")](https://rawgit.com/DMRobertson/thompsons_v/master/docs/examples/tree_pair/TreePair_render.svg "An element of F.")

Tree pair objects represent certain maps from an interval to itself. To plot the graph that this represents, run:

	>>> x.render_bijection(filename='example_plot.svg')

![TreePair("11000", "10100", "1 2 3") as a bijection](https://rawgit.com/DMRobertson/thompsons_v/master/docs/examples/tree_pair/TreePair_render_bijection.svg "The function that this element represents. It maps subintervals linearly: [0, 0.25] -> [0, 0.5]; [0.25, 0.5] -> [0.5, 0.75]; [0.5, 1] -> [0.75, 1].")

Tree Pairs can be multiplied to represent the composition of functions. They can also be inverted.
