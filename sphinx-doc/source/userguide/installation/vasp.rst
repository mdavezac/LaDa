.. _install_vasp_ug:

Setting up LaDa to call VASP
============================

There are only two variable specific to vasp calculations:
   
 - :py:data:`~lada.is_vasp_4` defines whether the installed vasp program is
   version 4.6 or 5.0 and higher. In practice, this determines which POSCAR
   format to use, and whether or not some input options are available.
 - :py:data:`~lada.vasp_program` defines the vasp executable. In general, it
   will be a string with path to the executable, although more powerful
   options are allowed.

.. note::

   LaDa should be :ref:`set up <install_mpi_ug>` properly to run mpi calculations.

.. warning::

   Please follow the links for their description. Questions regarding how to
   compile should be addressed to the relevant authorities.
