.. _install_mpi_ug:

Setting up MPI
==============


single call to mpi program
--------------------------


The command-line for mpi executables is stored in :py:data:`lada.mpirun_exe`::

  # openmpi and friends
  mpirun_exe = "mpirun -n {n} {placement} {program}"
  # Crays
  mpirun_exe = "aprun -n {n} {placement} {program}"

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
                                                     

.. _install_mmpi_ug:

multiple mpi programs in parallel
---------------------------------

It can be adavantageous to make simultaneous calls to different mpi programs or
the different calculations using the same mpi program. For instance, in a
genetic algorithm search, many different individuals need to be evaluated. It
would be best to do this in parallel and asynchronically. LaDa does provide for
this ability. However, setting it up can be difficult. The trick is to make
sure that LaDa knows which nodes it is allowed to run on each time it is
launched within a PBS or SLURM script.

LaDa figures out which nodes it can run on by launching at the start of the
calculation a small mpi program with all the allocated nodes, using 
directive defined by :py:data:`~lada.mpirun_exe`. This small program merely
calls socket.gethostname_ from each process. This way, it figures out the
nodes, and how many procs should be placed on each note. Later on, the
"{placement}" option in :py:data:`~lada.mpirun_exe` is used to specify where to
run each executable.

Lets figure out the different steps to setting this up. The easiest way to
figure out the value of a specific configuration variable is start an ipython
session and simply check what it is in the root :py:mod:`lada` package.

  1. Check that :py:data:`~lada.do_multiple_mpi_programs` is set to ''True''
  2. Check :py:data:`~lada.mpirun_exe` makes sense. Make sure that you
     have "{placement}" in this string (see step 4).
  3. Check :py:data:`~lada.figure_out_machines` and make sure that it will work
     with your installation. By default, it requires that boost.mpi_ and its
     boost.python_ bindings be installed. It should be fairly simple to change
     it to some other mpi python-bindings.
  4. If you are not using a standard mpich-like library (e.g. Cray crap), you
     will have to modify :py:func:`~lada.modify_global_comm` to suit your
     needs. Good luck with that machine.
  5. Check that :py:func:`~lada.placement` makes sense. This is a function
     which takes a :py:class:`~lada.process.mpi.Communicator` instance and returns a
     string which will be substituted for "{placement}" in
     :py:data:`~lada.mpirun_exe`. By default, it returns a path to a temporary
     `file listing the machines`_ to use. This should work on any supercomputer
     with a standard mpi installation. A communicator is a standard python
     dictionary with an additional
     :py:attr:`~lada.process.mpi.Communicator.machines` attribute. For standard
     mpi installations, you will only need to make sure that the correct option
     is used to specify a machine file (generally, '-machinefile').
  6. Check :py:data:`~lada.pbs_string` 

Often time, you will only need to go through steps 1 and 2, and install the
boost.mpi_ `python bindings <boost.python>`_. 

.. _boost.mpi: http://www.boost.org/doc/libs/1_49_0/doc/html/mpi.html
.. _boost.python: http://www.boost.org/doc/libs/1_49_0/doc/html/mpi/python.html
.. _format string: http://docs.python.org/library/st dtypes.html#str.format
.. _socket.gethostname: http://docs.python.org/library/socket.html#socket.gethostname
.. _file listing the machines: https://www.google.co.uk/search?q=mpi+machine+files
                                                     
                                                     
                                                     
                                                     
                                                     
