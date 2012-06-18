.. _lada-config:

Configuring LaDa
****************

Environment Variables
=====================

.. envvar:: LADA_CONFIG_DIR

   Environment variable specifying the path(s) to the configuration directories.

.. envvar:: LADA_DATA_DIRECTORY

   Environment variable specifying the path to the root of the data directories.

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
     
     .. note::

        Only taken into account at ipython start-up. It is ignored if LaDa is
        launched within python.

VASP 
----

  These variables are generally declared in config/vasp.py

  .. py:data:: is_vasp_4
     
     If it exists and is True, some vasp parameters will fail if used with
     vasp-5 only options. If it does not exist or is false, then these
     parameters are allowed. 

  .. py:data:: vasp_program

     Signifies which vasp executable to use. It can take the following values:

     - string: Should be the path to the vasp executable. It can be either
       a full path, or an executable within the envirnoment's $PATH
       variable.
     - callable: The callable is called with a
       :py:class:`~lada.vasp.functional.Vasp` instance as sole argument. It
       should return a string, as described above.  In other words, different
       vasp executables can be used depending on the parameters. 

       For instance, the following function chooses between a *normal* vasp and
       vasp compiled for perturbative spin-orbit calculations::

         def vasp_program(vasp):
           """ Signifies the vasp executable. 
           
               It is expected that two vasp executable exist, a *normal* vasp, and a one
               compiled for non-collinear calculations.
       
           """
           return "vasp-4.6-nc" if getattr(vasp, 'lsorbit', False) == True else "vasp-4.6"

.. _mpi-config:

MPI
---

  These variables are generally declared in config/mpi.py

  .. py:data:: mpirun_exe 
   
     Format string to launch mpi programs. It accepts as arguments 
     ``program``, ``n``, ``ppn`` as well as anything you want
     to throw at it:

     - ``program``: program (with commandline arguments) to launch
     - ``placement``: used to place an executable on specific nodes and processors
     - ``n``: number of processes to launch program
     - ``ppn``: number of processes per nodes
     
     In general, it takes the following form::

       mpirun_exe = "mpirun -n {n} {placement} {program}"

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

  .. py:data:: do_multiple_mpi_programs

     A boolean defining whether to attempt to figure out which machines lada
     can run on. This is only necessary if you will run different mpi
     executable simultaneously in the same PBS job.
  
  .. py:data:: figure_out_machines

     A string which defines a python program to get the hostnames of the
     machines LaDa can run on. This program must print to the standard output
     the names of the machines, one per line, and nothing else. Defaults to::

       figure_out_machines =  'from socket import gethostname\n'                      \
                              'from boost.mpi import world\n'                         \
                              'for i in xrange(world.size):\n'                        \
                              '  if i == world.rank: print gethostname(), i\n'        \
                              '  world.barrier()\n'                            

  .. py:function:: modify_global_comm(lada.process.mpi.Communicator)->None

     Called after figuring the hostnames of the nodes LaDa should run on. It is
     a callable tacking the global communicator as sole input. It should modify
     the callable such that :py:data:`~lada.placement` can make sense of a
     communicator and issue the correct placement configuration to the mpirun
     program. By default, this function does nothing.

  .. py:function:: placement(lada.process.mpi.Communicator)->str

     Callable which takes an :py:class:`~lada.mpi.Communicator` and returns a
     string which tells the mpirun program which nodes to run on. The string is
     substituted for "{placement}" in :py:data:`~lada.mpirun_exe`. In most
     cases (e.g. default), this means writing a machine file to disk and
     telling mpirun where it is with "-machinefile". 

Job-folder
----------

  .. py:data:: jobparams_readonly
   
     Whether instances of
     :py:class:`~lada.jobfolder.forwarding_dict.ForwardingDict` are read only
     by default. In practice, objects which use forwarding dictionaries
     generally dictate whether it should read-only or not, depending on what
     these objects do. This parameter should presently not have any effect.


  .. py:data:: jobparams_naked_end

     Whether mass collectors and manipulators, such as
     :py:class:`~lada.jobfolder.manipulator.JobParams` should return an object
     as is, rather than a
     :py:class:`~lada.jobfolder.forwarding_dict.ForwardingDict`, when it is the
     only item left. Practical when checking results in ipython, not so much
     when writing scripts.

  .. py:data:: jobparams_only_existing
    
     Whether, when setting parameters with
     :py:class:`~lada.jobfolder.manipulator.JobParams`, new attributes should
     be created for those items which do not possess that attribute, or whether
     :py:class:`~lada.jobfolder.manipulator.JobParams` should content itself
     with only modifying pre-existing attributes. Beware if set to True.
   
  .. py:data:: unix_re

      Whether mass collectors and manipulators, such as
      :py:class:`~lada.jobfolder.manipulator.JobParams`, accept regex as
      indices, or whether to use bash-like substitutions. The former is more
      powerfull, and the latter much simpler.

Compuational ressources and job submission
------------------------------------------

  .. py:data:: pbs_string

     String from which to create pbs_/slurm_ submission scripts. For instance,
     the following is for the slurm_ ressource manager::

        "#! /bin/bash/\n"                  \
        "#SBATCH --account={account}\n"    \
        "#SBATCH --time={walltime}\n"      \
        "#SBATCH -N {nnodes}\n"            \
        "#SBATCH -e {err}.%j\n"            \
        "#SBATCH -o {out}.%j\n"            \
        "#SBATCH -J {name}\n"              \
        "#SBATCH -D {directory}\n\n"       \
        "python {scriptcommand}\n"

     There are number of keywords which should appear:

        - walltime: defines how long the job should run. It will generally be
          provided when calling launch in ipython.
        - n: The number of processes to request from the resource manager.
        - nnodes: The number of nodes to request from the resource manager.
          Generally, it will be generated automatically from ``n`` and
          :py:data:`default_pbs`'s relevant information.
        - err: A file where to log errors from this job. This filename will be
          generated automatically.
        - out: A file where to log output from this job. This filename will be
          generated automatically.
        - name: The name of the job. Also generated automatically.
        - directory: The directory where the job will take place. Also
          generated automatically.
        - scriptcommand: You do want something to happen, right? Generated automatically.
        - account: Relevant to slurm_ only. Selected by user when launching job.

     Any number of parameters can be further provided, as long as they exist in
     :py:data:`default_pbs`.

  .. py:data:: default_pbs

     A dictionary which contains the parameters relevant to :py:data:`pbs_string`.
     Additionally, it should contain:
      
        - ppn: Number of processes per node. 

  .. py:data:: debug_queue
     
     How to select the debug queue. First part of the tuple is the keyword
     argument to modify when calling the pbs job, and the second is its value.

  .. py:function:: ipython_qstat

     An ipython magic function which returns all jobs submitted by the user.
     Once provided, it will be automatically imported into the ipython session
     by the lada extension, where is called ``qstat``.  This will change
     somewhat from one supercomputer to the next, depending on the type of
     ressource manager it uses. Here is what the function looks like for
     slurm_::

        def ipython_qstat(self, arg):
          """ squeue --user=`whoami` -o "%7i %.3C %3t  --   %50j" """
          from subprocess import Popen, PIPE
          from IPython.utils.text import SList
          from getpass import getuser
      
          # finds user name.
          whoami = getuser()
          squeue = Popen(["squeue", "--user=" + whoami, "-o", "\"%7i %.3C %3t    %j\""], stdout=PIPE)
          result = squeue.stdout.read().rstrip().split('\n')
          result = SList([u[1:-1] for u in result[1:]])
          return result.grep(str(arg[1:-1]))

     An this one is for the PBSpro_ ressource managers::

       def ipython_qstat(self, arg):
         """ Prints jobs of current user. """
         from subprocess import Popen, PIPE
         from IPython.utils.text import SList
         # get user jobs ids
         jobs   = SList(Popen(['qstat', '-f'], stdout=PIPE)                           \
                       .communicate()[0].split('\n'))
         names  = [ u[u.find('=')+1:].lstrip().rstrip()                               \
                    for u in jobs.grep('Job_Name') ]
         mpps   = [int(u[u.find('=')+1:]) for u in jobs.grep('Resource_List.ncpus')]
         states = [ u[u.find('=')+1:].lstrip().rstrip()                               \
                    for u in jobs.grep('job_state') ]
         ids    = [u[u.find(':')+1:].lstrip().rstrip() for u in jobs.grep('Job Id')]
         return SList([ "{0:>10} {1:>4} {2:>3} -- {3}".format(id, mpp, state, name)   \
                        for id, mpp, state, name in zip(ids, mpps, states, names)])

     It use IPython's SList_ to easily grep through output.

     Other/better snippets for other ressource managers are welcome.

  .. py:data:: queues

     List of strings defining the queues accessible to the users. They will be
     made available in :ref:`%lauch <ipython_launch_ug>`. It can be an empty
     tuple if "queues" are not relevant to the ressource manager.

  .. py:data:: accounts

     List of strings defining the accounts accessible to the users. They will
     be made available in :ref:`%lauch <ipython_launch_ug>`. It can be an empty
     tuple if "accounts" are not relevant to the ressource manager.

 .. _slurm: https://computing.llnl.gov/linux/slurm/
 .. _PBSpro: http://www.pbsworks.com/Product.aspx?id=1
 .. _SList: http://ipython.scipy.org/Wiki/Cookbook/StringListProcessing
 .. _pbs: http://www.mcs.anl.gov/research/projects/openpbs/
