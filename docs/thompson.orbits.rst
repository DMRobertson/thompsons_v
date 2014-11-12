Orbits
======

Let :math:`y` be a word, and let :math:`X` be an expansion of the standard basis :math:`\boldsymbol{x}=\{x_1, \dotsc, x_n\}`; finally let :math:`\phi` be an :class:`~thompson.automorphism.Automorphism`. We call the set :math:`\mathrm{orb}_\phi(y) = \{ y \phi^i\}_{i\in \mathbb Z}` the :math:`\phi`-orbit of :math:`y`. Higman demonstrated how we could classify these orbits in [Hig]_ (section 9). 

This module is responsible for two orbit-related data structures. Firstly, :class:`~thompson.orbits.OrbitType` is set up to describe all possible types of orbits according to Higman's scheme. Secondly, the :class:`~thompson.orbits.SolutionSet` describes solutions to the equation :math:`u \phi^i = v`.

The OrbitType structure
-----------------------

Let us refer to the intersection :math:`\mathrm{orb}_\phi (y) \cap X\langle A \rangle` as the :math:`X` *-component* of (the :math:`\phi`-orbit of) :math:`y`. These components come in five different types:
	
	1. Complete infinite.
		The component is the entire orbit :math:`\{ y \phi^i\}_{i\in \mathbb Z}` and each element of the orbit is different.
	2. Complete finite (or *periodic*).
		The component is the entire orbit, but eventually it repeats itself. Thus the component only contains a finite number of distinct words.
	3. Right semi-infinite.
		The forward orbit :math:`\{ y \phi^n\}_{n\in \mathbb N}` belongs to the component and does not repeat itself. However, the backward orbit :math:`\{ y \phi^{-n}\}_{n\in \mathbb N}` eventually leaves :math:`X\langle A\rangle`; thus the component only contains a finite amount of the backward orbit.
	4. Left semi-infinite.
		The backward orbit is contained in the component, and does not repeat itself; the forward orbit eventually leaves :math:`X\langle A\rangle`.
	5. Incomplete
		Only a finite part of the orbit :math:`y\phi^{-n}, \dotsc, y\phi^{-1}, y, y\phi, \dotsc, y\phi^m` belongs to the component.

Now suppose :math:`X` is a semi-normal basis for some automorphism. The orbits which intersect 

Orbits which are not incomplete are given a second, independant label.

	A.
		Type A components are those with characteristic :math:`(m, \epsilon)`, where :math:`\epsilon` is the empty word. These are precisely the periodic orbits.
	B.
		Type B components are those with characteristic :math:`(m, \Gamma)`, where :math:`\Gamma \neq \epsilon` is nontrivial. Type B components must be semi-infinite, **but not all semi-infinite components are of type B!** 
	C.
		Type C components are those which are not of type A or B---that is, they are the components which do not have a characteristic. If :math:`x` is of type C, then there exists :
		
		Complete infinite components are always of type C, but **some semi-infinite components are of type B too**.
	

.. autoclass:: thompson.orbits.OrbitType
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