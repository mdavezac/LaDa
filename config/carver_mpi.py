if "jobs" in globals()["pyladamodules"]:
  queues = "debug", "regular", "low"
  """ List of slurm or pbs queues allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """
