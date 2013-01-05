if "jobs" in globals()["pyladamodules"]:
  queues = "debug", "regular", "low", "premimum"
  """ List of slurm or pbs queues allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """
  resource_string = "mppwidth={nprocs}"
  """ Format string to specify computational resources. 
      
      The first argument is total number of processes, the second the number of
      nodes itself, the third the number of processes per node.
  """
  mpirun_exe = "aprun -n {nprocs} {program}"
  """ Command-line to launch external mpi programs. """
  cpus_per_node = 24
  """ Number of cpus per node. """

try_import_matplotlib = False
