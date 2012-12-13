Crystal structure à la CRYSTAL
==============================

There are presently three different input structure declarations:

  - :py:class:`~lada.dftcrystal.molecule.Molecule` defines a molecule, e.g. a
    0-d object.

  - :py:class:`~lada.dftcrystal.crystal.Crystal` defines a crystal, e.g. a
    3-d object.

  - :py:class:`~lada.dftcrystal.external.External` wraps around a
    :py:class:`~lada.crystal.cppwrappers.Structure` instance. It allows mixing
    both declarative and transformative structural definitions.

.. automodule:: lada.dftcrystal.molecule

.. autoclass:: Molecule
   :show-inheritance:
   :members: 
   :inherited-members:
   :exclude-members: raw, count, read_input, output_map, keyword, index

   Usage is similar to :py:class:`~lada.dftcrystal.crystal.Crystal`. Please
   look there for more information.

.. automodule:: lada.dftcrystal.crystal

.. autoclass:: Crystal
   :show-inheritance:
   :members: 
   :inherited-members:
   :exclude-members: raw, count, read_input, output_map, keyword, index

   .. attribute:: symmgroup

      Index or name of the space-group. Either form need make sense to
      CRYSTAL_. 

   .. attribute:: params

      List of crystal parameters. These are printed as are and in the same
      order directly to CRYSTAL_'s input.

   .. attribute:: atoms
       
      List of atomic sites in the initial structure. The items should be of
      type :py:class:`~lada.crystal.cppwrappers.Atom`. The easiest approach is
      to add them using :py:meth:`add_atom`.

.. automodule:: lada.dftcrystal.external

.. autoclass:: External
   :show-inheritance:
   :members: 
   :inherited-members:
   :exclude-members: raw, count, read_input, output_map, keyword, index
