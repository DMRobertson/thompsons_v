[Nathan Barker](https://www.dpmms.cam.ac.uk/~nb443/), [Andrew Duncan](http://www.mas.ncl.ac.uk/~najd2/) and [David Robertson](https://DMRobertson.github.io) are currently preparing a paper entitled *The power conjugacy problem in Higman-Thompson groups* which describes an algorithm to solve a certain equation in the groups named $G_{n,r}$.

This package aims to implement the algorithms described in the paper. To do so it provides tools for working in

- the algebra $V_{n,r}$, a certain set of words, and
- the automorphism group $G_{n,r} = \Aut V_{n,r}$.


Installation
------------

Either clone the repository with

	git clone https://github.com/DMRobertson/thompsons_v.git

or download [the ZIPped version](https://github.com/DMRobertson/thompsons_v/archive/master.zip) and extract it. For more detailed instructions, see [the documentation](http://thompsons-v.readthedocs.org/en/master/introduction.html#installation).

Documentation
-------------

Documentation is automatically built by [Read the Docs](http://thompsons-v.readthedocs.org/).

The HTML documentation is generated using [Sphinx](https://pypi.python.org/pypi/Sphinx). To build, ensure that Sphinx is installed (``pip install sphinx``), then navigate to [the docs folder](/docs) and run ``make html``. More interestingly, I'm using Sphinx's [doctest extension](http://sphinx-doc.org/ext/doctest.html) to run some simple automated tests. To use it, ensure Sphinx is installed and run ``make doctest`` from the docs directory.
