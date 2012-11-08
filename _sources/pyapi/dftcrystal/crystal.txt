Crystal structure à la CRYSTAL
==============================

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
