if "ladabase" in globals()["ladamodules"]:
  local_push_dir = "/projects/nrel/cid/database_tmp"
  """ Local directory where database stuff is pushed. """
if "jobs" in globals()["ladamodules"]:
  template_pbs = globals()["default_slurm"]
  """ Template pbs script to use. Depends on machine. """
  debug_queue = "queue", "inter"
  """ How to select the debug queue. 

      First part of the tuple is the keyword argument to modify when calling
      the pbs job, and the second is its value.
  """
  accounts = ["BES000"]
  """ List of slurm or pbs accounts allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """
  qsub_exe = "sbatch"
  """ Qsub executable. """
  resource_string = "-N {nnodes}"
  """ Format string to specify computational resources. 
      
      The first argument is total number of processes, the second the number of
      nodes itself, the third the number of processes per node.
  """
  mpirun_exe = "mpirun -np {nprocs} numa_wrapper -ppn={ppernode} {program}"
  """ Command-line to launch external mpi programs. """
  cpus_per_node = 8
  """ Number of processes per node. """
  lada_with_slurm = True
  """ If True use slurm as ressource manager, else use openpbs. """
