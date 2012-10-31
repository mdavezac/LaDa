mpirun_exe = "mpirun -n {n} {placement} {program}"
""" Command-line to launch external mpi programs. """
def machine_dependent_call_modifier(formatter=None, comm=None, env=None):
  """ Machine dependent modifications. 
  
      This is a fairly catch all place to put machine dependent stuff for mpi
      calls, including mpi placement.

      The formatter used to format the :py:data:`~lada.mpirun_exe` string is
      passed as the first argument. It can be modified *in-place* for machine
      dependent stuff, or for mpi placement. The latter case occurs only if
      ``comm`` has a non-empty ``machines`` attribute. In that case,
      :py:attr:`~lada.process.mpi.machines` is a dictionary mapping the
      hostnames to the number of procs on that host. Finally, an dictionary
      containing the environment variables can also be passed. It should be
      modified *in-place*.

      By default, the 'placement' value of the formatter is modified to reflect
      the nodefile of a specific mpi placement. This occurs only if
      mpi-placement is requested (eg `comm.machines` exists and is not empty).

      This function is called only from :py:function:`lada.launch_program`. If
      calls fail, it is a good idea to copy :py:function:`lada.launch_program`
      into your $HOME/.lada and debug it from there.

      :param dict formatter:
        Dictionary used in formatting the command line of
        :py:function:`~lada.launch`. It should be modified *in-place*.
      :param comm:
        Communicator used in this particular calculation. At this point in
        :py:function:`~lada.launch_program`, dictionary data from the
        communicator have been copied to the formatter. It is passed here in
        case its attributes :py:attr:`~lada.process.mpi.Communicator.machines`
        or the nodefile returned by
        :py:method:`~lada.process.mpi.Communicator.nodefile`
        is needed. However, the communicator itself should not be modified.
      :type comm: :py:class:`~lada.process.mpi.Communicator`
      :param dict env: 
        Dictionary of environment variables in which to run the call.

      :return: ignored
  """
  if len(getattr(comm, 'machines', [])) != 0:
    formatter['placement'] = "-machinefile {0}".format(comm.nodefile())

def modify_global_comm(communicator):
  """ Modifies global communicator so placement can be done correctly. 
  
      This function is called by :py:func:`create_global_comm`. It can be used
      to modify the global communicator to work better with a custom placement
      function.
  """
  pass

def launch_program( cmdl, comm=None, formatter=None, env=None, 
                    stdout=None, stderr=None, stdin=None, outdir=None ):
  """ Command used to launch a program.
  
      This function launches external programs for LaDa. It is included as a
      global so that it can be adapted to different computing environment. It
      also makes it easier to debug LaDa's mpi configuration when installing on
      a new machine.

      .. note::

        The number one configuration problem is an incorrect
        :py:data:`~lada.mpirun_exe`.

      .. note::
       
        The number two configuration problem is mpi-placement (eg how to launch
        two different mpi program simultaneously in one PBS/SLURM job). First
        read the manual for the mpi environment on the particular machine LaDa
        is installed on. Then adapt
        :py:function:`~lada.machine_dependent_call_modifier` by redeclaring it
        in $HOME/.lada.

      :param str cmld: 
        Command-line string. It will be formatted using ``formatter`` or
        ``comm`` if either are present. Otherwise, it should be exactly the
        (bash) command-prompt.
      :param comm: 
        Should contain everythin needed to launch an mpi call. 
        In practice, it is copied from :py:data:`~lada.default_comm` and
        modified for the purpose of a particular call (e.g. could use fewer
        than all available procs)
      :type comm: :py:class:`~lada.process.mpi.Communicator`
      :param dict formatter:
        Dictionary with which to format the communicator. If ``comm`` is
        present, then it will be updated with ``comm``'s input.
      :param dict env: 
        Dictionary containing the environment variables in which to do call.
      :param stdout:
        File object to which to hook-up the standard output. See Popen_.
      :param stderr:
        File object to which to hook-up the standard error. See Popen_.
      :param str outdir:
        Path to the working directory.

      .. _Popen:: http://docs.python.org/library/subprocess.html#popen-constructor
  """
  from shlex import split as shlex_split
  from subprocess import Popen
  from lada import machine_dependent_call_modifier
  from lada.misc import Changedir

  # make sure that the formatter contains stuff from the communicator, eg the
  # number of processes.
  if comm is not None and formatter is not None:
    formatter.update(comm)
  # Stuff that will depend on the supercomputer.
  machine_dependent_call_modifier(formatter, comm, env)

  # if a formatter exists, then use it on the cmdl string.
  if formatter is not None: cmdl = cmdl.format(**formatter)
  # otherwise, if comm is not None, use that.
  elif comm is not None: cmdl = cmdl.format(**comm)

  # Split command from string to list
  cmdl = shlex_split(cmdl)

  # makes sure the directory exists:
  if outdir is not None:
    with Changedir(outdir) as cwd: pass
  # finally, start process.
  return Popen( cmdl, stdout=stdout, stderr=stderr, stdin=stdin, cwd=outdir,
                env=env )



default_comm = {'n': 2, 'ppn': 4, 'placement': ''}
""" Default communication dictionary. 

    should contain all key-value pairs used in :py:data:`mpirun_exe`.  In a
    script which manages mpi processes, it is also the global communicator. In
    other words, it is the one which at the start of the application is given
    knowledge of the machines (via :py:func:`~lada.create_global_comm`). Other
    communicators will have to acquire machines from this one. In that case, it
    is likely that 'n' is modified.
"""

# pbs/slurm related stuff.
queues = ()
""" List of slurm or pbs queues allowed for use. 

    This is used by ipython's %launch magic function. 
    It is not required for slurm systems. 
    If empty, then %launch will not have a queue option.
"""
accounts = ['BES000']
""" List of slurm or pbs accounts allowed for use. 

    This is used by ipython's %launch magic function. 
    It is not required for slurm systems. 
    If empty, then %launch will not have a queue option.
"""

debug_queue = "queue", "debug"
""" How to select the debug queue. 

    First part of the tuple is the keyword argument to modify when calling
    the pbs job, and the second is its value.
"""
qsub_exe = "sbatch"
""" Qsub/sbatch executable. """
qsub_array_exe = None
""" Qsub for job arrays.

    If not None, if should be a tuple consisting of the command to launch job
    arrays and the name of the environment variable holding the job index. 

    >>> qsub_array_exe = 'qsub -J 1-{nbjobs}', '$PBS_ARRAY_INDEX'

    The format ``{array}`` will receive the arrays to launch.
"""
qdel_exe = 'scancel'
""" Qdel/scancel executable. """

default_pbs = { 'account': accounts[0], 'walltime': "06:00:00", 'nnodes': 1,
                'ppn': 1, 'header': "", 'footer': "" }
""" Defaults parameters filling the pbs script. """
pbs_string =  "#! /bin/bash/\n"                                                \
              "#SBATCH --account={account}\n"                                  \
              "#SBATCH --time={walltime}\n"                                    \
              "#SBATCH -N={nnodes}\n"                                          \
              "#SBATCH -e={err}\n"                                             \
              "#SBATCH -o={out}\n"                                             \
              "#SBATCH -J={name}\n"                                            \
              "#SBATCH -D={directory}\n\n"                                     \
              "{header}\n"                                                     \
              "python {scriptcommand}\n"                                       \
              "{footer}\n"
""" Default pbs/slurm script. """

do_multiple_mpi_programs = True
""" Whether to get address of host machines at start of calculation. """

figure_out_machines =  'from socket import gethostname\n'                      \
                       'from boost.mpi import gather, world\n'                 \
                       'hostname = gethostname()\n'                            \
                       'results = gather(world, hostname, 0)\n'                \
                       'if world.rank == 0:\n'                                 \
                       '  for hostname in results:\n'                          \
                       '    print "LADA MACHINE HOSTNAME:", hostname\n'        \
                       'world.barrier()\n'
""" Figures out machine hostnames for a particular job.

    Can be any programs which outputs each hostname (once per processor),
    preceded by the string "LADA MACHINE HOSTNAME:"
"""
