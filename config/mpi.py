mpirun_exe = "mpirun -n {n} {placement} {program}"
""" Command-line to launch external mpi programs. """
def placement(communicator=None):
  """ Placement string for MPI processes. 
  
      Should return an empty string if the communicator is None.
      This version works for openmpi.
  """
  if communicator is None: return ""
  if len(getattr(communicator, 'machines', {})) == 0: return ""
  return "-machinefile {0}".format(communicator.nodefile())

def modify_global_comm(communicator):
  """ Modifies global communicator so placement can be done correctly. 
  
      This function is called by :py:func:`create_global_comm`. It can be used
      to modify the global communicator to work better with a custom placement
      function.
  """
  pass


default_comm = {'n': 2, 'ppn': 4}
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

default_pbs = {'account': accounts[0], 'walltime': "06:00:00", 'nnodes': 1 }
""" Defaults parameters filling the pbs script. """
pbs_string =  "#! /bin/bash/\n"\
              "#SBATCH --account={account}\n"\
              "#SBATCH --time={walltime}\n"\
              "#SBATCH -N={nnodes}\n"\
              "#SBATCH -e={err}\n"\
              "#SBATCH -o={out}\n"\
              "#SBATCH -J={name}\n"\
              "#SBATCH -D={directory}\n\n"\
              "python {scriptcommand}\n"
""" Default pbs/slurm script. """

do_multiple_mpi_programs = True
""" Whether to get address of host machines at start of calculation. """

figure_out_machines =  'from socket import gethostname\n'                      \
                       'from boost.mpi import world\n'                         \
                       'for i in xrange(world.size):\n'                        \
                       '  if i == world.rank: print gethostname(), i\n'        \
                       '  world.barrier()\n'                            

