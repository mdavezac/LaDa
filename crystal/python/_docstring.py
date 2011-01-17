""" Contains basic data type and methods for crystal structure and lattices.

    Example: Creating a silicon nanowire inside a germanium structure
    =================================================================

    The first step will be to create the Si/Ge zinc-blende lattice. The lattice
    object will allow us to create any kind of supercell with a single line of
    code. From there, we need only loop over all atoms in our super-structure
    and assign the type ("Si" or "Ge") depending on whether it is inside or
    outside the nanowire.

    Creating a Si/Ge zinc-blende lattice
    ------------------------------------

    There are two main types of objects in `lada.crystal`: `Lattice`, wich
    defines the unit-cell of a crystal structure, and `Structure` which defines
    any kind of crystal structure, including super-cells. The advantage of
    difining the lattice as a separate object is that any super-structure can
    then be created very easily.

    The zinc-blende lattice has already been defined by default. It can be
    accessed simply as:

    >>> from lada.crystal import binary
    >>> zincblende = binary.zinc_blende()
    >>> zincblende
    # Lattice definition.
    from lada.crystal._crystal import Lattice
    lattice = Lattice()
    lattice.name   = 'Zinc-Blende'
    lattice.scale  = 1.0
    lattice.set_cell = (0.0, 0.5, 0.5),\\
                       (0.5, 0.0, 0.5),\\
                       (0.5, 0.5, 0.0)
    lattice.add_sites = [(0.00, 0.00, 0.00), 'A'],\\
                        [(0.25, 0.25, 0.25), 'B']
    # End of lattice definition

    The first line above import a module containing all sorts of binary
    lattices. The second line calls a function which returns a zinc-blende
    lattice and stores it in the variable ``zinc_blende`` for further
    processing. This last object can be printed (third line): it outputs python
    code which defines the lattice itself. This code is fairly
    self-explanatory, and can be copied/pasted and modified to suit your needs
    if a different kind of lattice is needed. 
    
    As everywhere in `lada`, atomic-positions are given in cartesian
    coordinates. The units are defined by both ``zincblende.scale`` and the
    cell and lattice positions. The scale is equivalent to the first line of a
    POSCAR file. Eventually, all units are in Angstrom.

    As can be seen, the atomic types are "A" and "B". We would rather have that
    each site can be occupied by either "Si" or "Ge":

    >>> for site in zincblende.sites:
    >>>   site.type = ["Si", "Ge"]
    >>> zincblende
    # Lattice definition.
    from lada.crystal._crystal import Lattice
    lattice = Lattice()
    lattice.name   = 'Zinc-Blende'
    lattice.scale  = 1.0
    lattice.set_cell = (0.0, 0.5, 0.5),\\
                       (0.5, 0.0, 0.5),\\
                       (0.5, 0.5, 0.0)
    lattice.add_sites = [(0.00, 0.00, 0.00), ('Si', 'Ge')],\\
                        [(0.25, 0.25, 0.25), ('Si', 'Ge')]
    # End of lattice definition.
    
    The first line above loops over all sites in the lattice. The second line
    the modifies the type of occupations which can occur in the lattice.

    Creating a superstructure
    -------------------------

    At this point, we have all the ingredients to create a superstructure which
    will hold a Silicon nanowire:

    >>> structure = zincblende.to_structure([[5,0,0], [0,2,0], [0, 0,2]]) 
    >>> structure
    # Structure definition.
    from lada.crystal._crystal import Structure
    structure = Structure()
    structure.name   = ''
    structure.scale  = 1.0
    structure.energy = 0.0
    structure.weight = 1.0
    structure.set_cell = (5.0, 0.0, 0.0),\\
                         (0.0, 2.0, 0.0),\\
                         (0.0, 0.0, 2.0)
    structure.add_atoms = [(0.00, 0.00, 0.00), 'Si', 0],\\
                          [(0.25, 0.25, 0.25), 'Si', 1],\\
                          ...
                          [(4.50, 1.50, 0.00), 'Si', 0],\\
                          [(4.75, 1.75, 0.25), 'Si', 1]
    # End of structure definition.

    The first line above creates a super-structure with an extension of 5
    conventional cells in the (100) direction, and extensions of 2 conventional
    cells in (010) and (001) directions.

    The atoms exist in the super-structure. However, the occupations all
    default to "Si" (since it is the first possible occupation of the lattice
    sites, as defined above).

    Creating the nanowire within the superstructure
    -----------------------------------------------

    The first step to creating a nanowire is making sure we won't have any
    problems with periodic images (They may, or may not be a this point,
    depending how the zinc-blende lattice was defined). We can move all atoms
    in the superstructure such that they are contained within its supercell. To
    do this, we go to fractional coordinates, then move the fractional
    coordinates to the interval [0, 1[, and finally move back to cartesian
    coordinates.

    >>> # some numpy functions we will need.
    >>> from numpy.linalg import inv
    >>> from numpy import dot, floor, all
    >>> # compute inverse of superstructure cell.
    >>> invcell = inv(structure.cell)
    >>> # Loop over all atoms in superstructure.
    >>> for atom in structure.atoms:
    >>>   # Compute fractional coordinates.
    >>>   fractional = dot(invcell, atom.pos)
    >>>   # Move fractional coordinates to [0, 1[ interval, taking into account numerical noise.
    >>>   fractional -= floor(fractional + 1e-12)
    >>>   # Make sure that we did as advertised. Always a good idea.
    >>>   assert all(fractional >= 0e0) and all(fractional < 1e0)
    >>>   # Move from fractional back to cartesian coordinates.
    >>>   atom.pos = dot(structure.cell, fractional)


    We can now change the occupation such that atoms within the nanowire a "Si", and 
    "Ge" outside. To do this, we simply compute the radial coordinates of the
    atom, where (1,0,0) is the growth direction.

    >>> # First, we define the growth direction.
    >>> from numpy import array
    >>> from numpy.linalg import norm
    >>> growth = array([1,0,0], dtype='float64') 
    >>> growth = growth / norm(growth) # make sure the vector is normalized
    >>> radius = 2.5 # just a random radius... in angstrom.
    >>> # loop over all atoms.
    >>> for atom in structure.atoms:
    >>>   # compute radial vector (eg remove growth component) of the atom.
    >>>   radial_vector = atom.pos - dot(growth, atom.pos)
    >>>   # compute radial norm
    >>>   radial_norm = norm(radial_vector)
    >>>   # assign type depending on whether atom is within radius or not.
    >>>   atom.type = "Si" if radial_norm < radius else "Ge"
    

    Coding Details
    ==============

    C++ bindings are located in the private `_crystal` bindics. A fair number of
    enhancements are added directly within the python code in __init__.py. In
    practice all public interfaces to C++ bindings should be available directly
    in the `crystal` module.
"""

