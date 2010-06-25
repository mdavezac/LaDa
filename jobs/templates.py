""" Templates for job submissal.

 
"""

def default_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = "regular", name = None, \
                 pyvirt = "vasp-4.6.32", pbspools = 1, npbs = 0, procpools = 1, \
                 pyscript = None, pickle = "job_pickle", outdir = None, **kwargs ):
  """ Creates default pbs-script. Does not launch. 

      @param file: File object to which to write. 
      @param walltime: must be a string in "hh:mm:ss" format, or anything
         acceptable to the PBS implementation.
      @param mppwidth: Number of processes (not processors) to use.
      @param queue: Queue to use
      @param name: Name of the job.
      @param pyvirt: Virtual python environment. Does nothing if None.
      @param pbspools: Number of pbs pools.
      @param npbs: Which pool is this script.
      @param procpools: Number of mpi pools (parallelization over processes).
      @param pyscript: Python script to launch. Default: L{runme.py}
      @param pickle: Fileame of job-tree pickle.
      @param kwargs: Ignored by this script.
  """
  from os.path import exists

  if name == None: name = "automatic" if pbspools <= 1 else "automatic_%i" % (npbs)

  file.write("#! /bin/bash\n")
  file.write("#PBS -l walltime=%s,mppwidth=%i\n" % (walltime, mppwidth) )
  file.write("#PBS -q " + queue + "\n")
  file.write("#PBS -e err\n")
  file.write("#PBS -o out\n")
  if name != None: file.write("#PBS -N %s \n\n" % (name))

  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd " + outdir + "\n")
  if pyvirt != None: file.write("workon %s \n" % (pyvirt) )

  file.write("mpirun -n %i python %s " % (mppwidth, pyscript) )
  if pbspools > 1: file.write("--pbspools %i --npbs %i" % (pbspools, npbs) )
  if procpools > 1: file.write("--procpools %i" % (procpools) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --%s" % (key))
    else:             file.write(" --%s %s" % (key, value))
  file.write(" " + pickle + "\n")
