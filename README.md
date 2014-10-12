[Nathan Barker](https://www.dpmms.cam.ac.uk/~nb443/) and [Andrew Duncan](http://www.mas.ncl.ac.uk/~najd2/) are currently preparing a paper entitled *The power conjugacy problem in Thompson's group V*, which describes an algorithm to solve the power conjugacy problem in the Higman-Thompson groups $G_{n,r}$.

This is a package for Python 3.3+ which aims to implement that algorithm. To do so it provides tools for working with the algebra $V_{n,r}$.

Installation
------------

Either clone the repository with

	git clone https://github.com/DMRobertson/thompsons_v.git

or download [the ZIPped version](https://github.com/DMRobertson/thompsons_v/archive/master.zip) and extract it.

Documentation
-------------

The HTML documentation is generated using [Sphinx](https://pypi.python.org/pypi/Sphinx). To build, ensure that Sphinx is installed (``pip install sphinx``), then navigate to [the docs folder](/docs) and run ``make html``. More interestingly, I'm using Sphinx's [doctest extension](http://sphinx-doc.org/ext/doctest.html) to run some simple automated tests. To use it, ensure Sphinx is installed and run ``make doctest`` from the docs directory.
