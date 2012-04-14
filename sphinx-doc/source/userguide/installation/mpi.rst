.. _install_mpi_ug:

Setting up MPI calls
====================


The command-line for mpi executables is stored in :py:data:`lada.mpirun_exe`::

  # openmpi and friends
  mpirun_exe = "mpirun -n {n} {program}"
  # Crays
  mpirun_exe = "aprun -n {n} {program}"

This is a `format string`_. It can take any number of arguments. However, two
are required: "n", which is the number of processes, and "program", which is
the commandline for the program to launch. The latter will be manufactured by
LaDa internally. It is a placeholder at this point. The other reseverved
keyword is "ppn", the number of processes per node. It should only be used for
that purpose. 

The keywords in :py:data:`lada.mpirun_exe` should be defined in
:py:data:`lada.default_comm`. This is a dictionary which holds default values
for the different keywords. The dictionary may hold more keywords than are
present in :py:data:`lada.mpirun_exe`. The converse is not true.
                                                     
.. _format string: http://docs.python.org/library/st dtypes.html#str.format
                                                     
                                                     
                                                     
                                                     
                                                     
