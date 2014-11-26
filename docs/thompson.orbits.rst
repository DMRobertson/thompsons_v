Orbits
======

Let :math:`y` be a word, and let :math:`X` be an expansion of the standard basis :math:`\boldsymbol{x}=\{x_1, \dotsc, x_n\}`; finally let :math:`\phi` be an :class:`~thompson.automorphism.Automorphism`. We call the set :math:`\mathrm{orb}_\phi(y) = \{ y \phi^i\}_{i\in \mathbb Z}` the :math:`\phi`-orbit of :math:`y`. Let us refer to the intersection :math:`\mathrm{orb}_\phi (y) \cap X\langle A \rangle` as the :math:`X` *-component* of (the :math:`\phi`-orbit of) :math:`y`. Higman demonstrated how we could classify these components in [Hig]_ (section 9). 

This module is responsible for two orbit-related data structures. Firstly, :class:`~thompson.orbits.OrbitType` is set up to describe all possible types of orbits according to Higman's scheme. Secondly, the :class:`~thompson.orbits.SolutionSet` describes solutions to the equation :math:`u \phi^i = v`.

Types of orbits and components
------------------------------

These components come in five different types:
	
	1. Complete infinite components.
		The component is the entire orbit :math:`\{ y \phi^i\}_{i\in \mathbb Z}` and each element of the orbit is different.
	2. Complete finite components.
		The component is the entire orbit, which eventually it repeats itself. Thus the component only contains a finite number of distinct words.
	3. Right semi-infinite components.
		The forward orbit :math:`\{ y \phi^n\}_{n\in \mathbb N}` belongs to the component and does not repeat itself. However, the backward orbit :math:`\{ y \phi^{-n}\}_{n\in \mathbb N}` eventually leaves :math:`X\langle A\rangle`; thus the component only contains a finite amount of the backward orbit.
	4. Left semi-infinite components.
		The backward orbit is contained in the component, and does not repeat itself; the forward orbit eventually leaves :math:`X\langle A\rangle`.
	5. Incomplete finite components.
		A finite part of the orbit
		
		.. math:: y\phi^{-n}, \dotsc, y\phi^{-1}, y, y\phi, \dotsc, y\phi^m
		
		for which both :math:`y\phi^{-n-1}` and :math:`y\phi^{m+1}` are not in :math:`X\langle A\rangle`.

There are six different ways we can form orbits using these components:
	
	1. Complete infinite (or *doubly infinite*) orbits.
		The entire orbit is a complete infinite component (type 1).
	2. Complete finite (or *periodic*) orbits.
		The entire orbit is a component finite component (type 2).
	3. Incomplete finite (or *bad*) orbits.
		The orbit does not contain any complete or semi-infinite components. Thus the only components which may appear are of type 5. Note that these orbits need not have any components at all!
	4. Incomplete infinite orbits.
		The orbit contains at least one semi-infinite component (types 3, 4). The orbit may also contain incomplete finite components (type 5). We give names to the different cases:
		
		a. Right semi-infinite orbits.
			One component of type 3, none of type 4.
		b. Left semi-infinite orbits.
			One component of type 4, none of type 3.
		c. *Pond* orbits.
			One component of type 3 and one of type 4. These components must be separated by a finite list of words in :math:`V_{n,r}\setminus X\langle A\rangle`; we think of this list as being a *pond*.

If that wasn't confusing enough, we have a another way to classify those orbits which are not incomplete finite.

	A.
		Type A components are those with characteristic :math:`(m, \epsilon)`, where :math:`\epsilon` is the empty word. These are precisely the periodic orbits.
	B.
		Type B components are those with characteristic :math:`(m, \Gamma)`, where :math:`\Gamma \neq \epsilon` is nontrivial. 
	C.
		Type C components are those which are not of type A or B---that is, they are the components which do not have a characteristic. If :math:`x` belongs to a type C orbit, then :math:`x\phi^\ell = z\Delta` for some :math:`z` of type B. That is, type C orbits are those 'below' type B orbits.
	
	.. note:: 
		
		Type B components must be semi-infinite, **but not all semi-infinite components are of type B**. Complete infinite components are always of type C, but **some semi-infinite components are of type C too**.
	

.. autoclass:: thompson.orbits.ComponentType
	:members:
	:undoc-members:

The SolutionSet structure
-------------------------

Solutions to the equation :math:`u\psi^m = v` come in specific instances. If there are no solutions, the solution set is :math:`\emptyset`. Otherwise :math:`u` and :math:`v` are distinct elements which share an orbit. If this orbit is periodic, then the solution set is :math:`m + p\mathbb{Z}`, where :math:`p` is the period of the orbit. Otherwise the solution set is a single integer :math:`m`.

Internally we represent this as a pair *(base, increment)* of integers. The first element *base* is a solution :math:`m` if it exists; otherwise ``None``. The second element *increment* is the increment *p* between solutions (which occurs for a periodic orbit only).

.. autoclass:: thompson.orbits.SolutionSet
	:members:
	:undoc-members:
	:exclude-members: singleton, empty_set, the_integers

Helper functions
----------------

.. automodule:: thompson.orbits
    :members:
    :undoc-members:
    :exclude-members: OrbitType, SolutionSet

Next steps
----------

With tools to describe orbits in hand, we can dive straight into the Automorphisms module. In particular, we need to :meth:`know how orbits work <thompson.automorphism.Automorphism.orbit_type>` to

- determine :meth:`quasinormal bases <thompson.automorphism.Automorphism.quasinormal_basis>`
- perform the :meth:`orbit sharing test <thompson.automorphism.Automorphism.share_orbit>`
- perform the :meth:`conjugacy test <thompson.automorphism.Automorphism.test_conjugate_to>`