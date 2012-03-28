.. _lada-config:

Configuring lada
****************

Environment Variables
=====================

.. envvar:: LADA_CONFIG_DIR

   Environment variable specifying the path(s) to the configuration directories.

Configuration variables
=======================

Configuration variables exist in the :py:mod:`lada` module itself. However,
they can be added within separate files. Which files will depend upon the user.

   - Files located in the config sub-directory where lada is installed
   - Files located in one of the directories specified by :envvar:`LADA_CONFIG_DIR`
   - In the user configuration file ~/.lada

Each file is executed and whatever is declared within is placed directly at the
root of the :py:mod:`lada` package. The files are read in that order. Within a
given directory, files are read alphabetically. Later files can override
previous files.

.. currentmodule:: lada

.. _vasp-config:

General
-------

  .. py:data:: verbose_representation
    
     Whether functionals should be represented/printed verbosely, e.g. each and
     every attribute, or whether attributes which have not changed from the
     default should be stripped. The former is safer since it should defaults
     may change over time, and the representation can become inaccurate.
     Defaults to False.
  
  .. py:data:: ipython_verbose_representation
    
     When in ipython and if not None, then changes
     :py:data:`verbose_representation` to this value. Makes it a bit easier on
     the eyes in ipython, while keeping things accurate during actual
     calculations. Ignored if None. Defaults to False. 
     
     .. note:: Only taken into account at ipython start-up.

VASP 
----

  These variables are generally declared in config/vasp.py

  .. py:data:: is_vasp_4
     
     If it exists and is True, some vasp parameters will fail if used with
     vasp-5 only options. If it does not exist or is false, then these
     parameters are allowed. 

  .. py:data:: vasp_program

     Path to the vasp executable itself.

.. _mpi-config:

MPI
---

  These variables are generally declared in config/mpi.py

  .. py:data:: mpirun_exe 
   
     Format string to launch mpi programs. It accepts as arguments 
     ``program``, ``commandline``, ``n``, ``ppn`` as well as anything you want
     to throw at it:

     - ``program``: program to launch
     - ``cmdline``: command line arguments to lauch it with
     - ``n``: number of processes to launch program
     - ``ppn``: number of processes per nodes
     
     In general, it takes the following form:

     >>> mpirun_exe = "mpirun -n {n} {program} {cmdline}

     The actual commandline is executed by :py:func:`execute_program
     <lada.misc.execute_program>`. The latter executes via Popen_ a
     commandline obtained through the format_ method of a python string. The
     arguments to format are those mentioned above as well as anything passed
     on to :py:func:`execute_program <lada.misc.execute_program>`.

     .. _Popen: http://docs.python.org/library/subprocess.html#subprocess.Popen
     .. _format: http://docs.python.org/library/stdtypes.html#str.format

  .. py:data:: default_comm 

     An dictionary with ``n`` and ``ppn``, as well as any other variable to be
     used in conjunction with :py:data:`mpirun_exe`.

