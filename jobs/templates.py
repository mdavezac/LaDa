""" Templates for job submissal.

 
"""

def default_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = "regular", name = None, \
                 pyvirt = "vasp.4.6.32", pools = 1,  pyscript = None, pickle = "job_pickle", \
                 outdir = None, python_path = None):
  """ Creates default pbs-script. Does not launch. 

      @param file: File object to which to write. 
      @param walltime: must be a string in "hh:mm:ss" format, or anything
         acceptable to the PBS implementation.
      @param mppwidth: Number of processes (not processors) to use.
      @param queue: Queue to use
      @param name: Name of the job.
      @param pyvirt: Virtual python environment. Does nothing if None.
      @param pools: Number of mpi pools (parallelization over processes).
      @param pyscript: Python script to launch. Default: L{runme.py}
      @param pythonpath: Additional path to search for modules.
      @param pickle: Fileame of job-tree pickle.
  """
  from os.path import exists

  file.write("#! /bin/bash\n")
  file.write("#PBS -l walltime=%s,mppwidth=%i\n" % (walltime, mppwidth) )
  file.write("#PBS -q " + queue + "\n")
  file.write("#PBS -e err\n")
  file.write("#PBS -o out\n")
  if name != None: file.write("#PBS -N %s \n\n" % (name))

  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd " + outdir + "\n")
  if pyvirt != None: file.write("workon %s \n" % (pyvirt) )

  file.write("mpirun -n %i python %s --pools %i " % (mppwidth, pyscript, pools) )
  if python_path != None: file.write(" --ppath " + python_path)
  file.write(" " + pickle + "\n")

