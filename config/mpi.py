mpirun_exe = "mpirun -n {n} {program}"
""" Command-line to launch external mpi programs. """

default_comm = {'n': 1, 'ppn': 8}
""" Default communication directory. 


    should contain all key-value pairs used in :py:data:`mpirun_exe`.
"""

if "jobs" in globals()["ladamodules"]:
  default_walltime = "06:00:00"
  """ Default walltime when launching jobs. """
  lada_with_slurm = False
  """ If True use slurm as ressource manager, else use openpbs. """
  queues = ()
  """ List of slurm or pbs queues allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """
  accounts = []
  """ List of slurm or pbs accounts allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """

  template_pbs = globals()["default_pbs"]
  """ Template pbs script to use. Depends on machine. """

  debug_queue = "queue", "debug"
  """ How to select the debug queue. 

      First part of the tuple is the keyword argument to modify when calling
      the pbs job, and the second is its value.
  """
  qsub_exe = "qsub"
  """ Qsub executable. """
  resource_string = "nodes={1}:ppn={2}" 
  """ Format string to specify computational resources. 
      
      The first argument is total number of processes, the second the number of
      nodes itself, the third the number of processes per node.
  """

  cpus_per_node = 4 #globals()["cpus_per_node"]()
  """ Number of processes per node. """
