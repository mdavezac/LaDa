""" Templates for job submissal.

 
"""

def default_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = "regular", name = None, \
                 pyvirt = None, pyscript = None, pickle = "job_pickle", \
                 outdir = None, **kwargs):
  """ Creates default pbs-script. Does not launch. 

      @param file: File object to which to write. 
      @param walltime: must be a string in "hh:mm:ss" format, or anything
         acceptable to the PBS implementation.
      @param mppwidth: Number of processes (not processors) to use.
      @param queue: Queue to use
      @param name: Name of the job.
      @param pyvirt: Virtual python environment. If None, checks for
        $VIRTUAL_ENV environment variable and uses that if it exists. Otherwise
        does nothing.
      @param pyscript: Python script to launch. Default: L{runme.py}
      @param pythonpath: Additional path to search for modules.
      @param pickle: Fileame of job-tree pickle.
  """
  from os import environ
  from os.path import exists, split as splitpath

  if pyvirt == None and "VIRTUAL_ENV" in environ:
    pyvirt = join(join(environ["VIRTUAL_ENV"], "bin"), "activate")


  file.write("#! /bin/bash\n")
  file.write("#PBS -l walltime=%s,mppwidth=%i\n" % (walltime, mppwidth) )
  file.write("#PBS -q " + queue + "\n")
  file.write("#PBS -e err.$PBS_JOBID\n")
  file.write("#PBS -o out.$PBS_JOBID\n")
  if name != None: file.write("#PBS -N %s \n\n" % (name))

  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd " + outdir + "\n")
# if pyvirt != None: file.write("\nsource %s \n" % (pyvirt) )

  # aprun on Fucking Crays. mpirun everywhere else.
  file.write("aprun " if "NERSC_HOST" in environ else "mpirun ")
  file.write("-n %i python %s " % (mppwidth, pyscript) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --%s" % (key))
    else:             file.write(" --%s %s" % (key, value))
  file.write(" " + pickle + "\n")

def default_slurm( file, walltime = "05:45:00", mppwidth = 8, ppernode=8, queue = None,\
                   name = None, pyvirt = None, pyscript = None,  pickle = "job_pickle",\
                   outdir = None, **kwargs):
  """ Creates default slurm-script. Does not launch. 

      @param file: File object to which to write. 
      @param walltime: must be a string in "hh:mm:ss" format, or anything
         acceptable to the PBS implementation.
      @param mppwidth: Number of processes (not processors) to use.
      @param ppernode: Number of processes per node.
      @param queue: Queue to use. If none will let slurm use default queue.
      @param name: Name of the job.
      @param pyvirt: Virtual python environment. If None, checks for
        $VIRTUAL_ENV environment variable and uses that if it exists. Otherwise
        does nothing.
      @param pyscript: Python script to launch. Default: L{runme.py}
      @param pythonpath: Additional path to search for modules.
      @param pickle: Fileame of job-tree pickle.
  """
  from os import environ
  from os.path import exists, split as splitpath, join

  if pyvirt == None and "VIRTUAL_ENV" in environ:
    pyvirt = join(join(environ["VIRTUAL_ENV"], "bin"), "activate")

  nnodes = mppwidth / ppernode + (0 if mppwidth % ppernode == 0 else 1)

  file.write("#! /bin/bash\n")
  file.write("#SBATCH --time=" + walltime + "\n") 
  file.write("#SBATCH -N " + str(nnodes) + "\n") 
  file.write("#SBATCH -n " + str(mppwidth) + "\n") 
  if queue != None: file.write("#SBATCH -p %s\n" % (queue))
  pbsdir = splitpath(file.name)[0]
  if name != None:
    file.write("#SBATCH -e \"%s/err.%s.%%j\"\n" % (pbsdir, name))
    file.write("#SBATCH -o \"%s/out.%s.%%j\"\n" % (pbsdir, name))
  else:
    file.write("#SBATCH -e \"%s/err.%j\"\n" % (pbsdir))
    file.write("#SBATCH -o \"%s/out.%j\"\n" % (pbsdir))
  if name != None: file.write("#SBATCH -J \"%s\" \n" % (name))
  if outdir != None: file.write("#SBATCH -D %s " % (outdir))

  if pyvirt != None: file.write("\nsource %s \n" % (pyvirt) )

  # aprun on Fucking Crays. mpirun everywhere else.
  file.write( "mpirun -np %i numa_wrapper -ppn=%i python %s "\
              % (mppwidth, ppernode, pyscript) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --%s" % (key))
    else:             file.write(" --%s %s" % (key, value))
  file.write(" " + pickle + "\n")

