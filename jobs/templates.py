""" Templates for job submissal. """
__docformat__ = "restructuredtext en"

def default_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = None, name = None, \
                 pyscript = None, pickle = "job_pickle", outdir = None, **kwargs):
  """ Creates default pbs-script. Does not launch. 

      :Parameters:
        file
          File object to which to write. 
        walltime
         must be a string in "hh:mm:ss" format, or anything
         acceptable to the PBS implementation.
        mppwidth
          Number of processes (not processors) to use.
        queue
          Queue to use
        name
          Name of the job.
        pyscript
          Python script to launch. Default: L{runme.py}
        pickle
          Fileame of job-tree pickle.
  """
  from os import environ
  from os.path import exists, join, dirname, abspath

  pbsdir = abspath(dirname(file.name))
  file.write("#! /bin/bash\n#PBS -l walltime={0},mppwidth={1}\n".format(walltime, mppwidth))
  if name != None: 
    file.write("#PBS -N {1}\n\n"\
               "#PBS -e {0}/err.{1}.$PBS_JOBID\n"\
               "#PBS -o {0}/out.{1}.$PBS_JOBID\n".format(pbsdir, name))
  else:
    file.write("#PBS -e {0}/err.$PBS_JOBID\n"\
               "#PBS -o {0}/out.$PBS_JOBID\n".format(pbsdir))
  if queue != None: file.write("#PBS -q {0} \n".format(queue))
  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd {0}\n".format(outdir))
  for key, value in environ.items():
    file.write("export {0}={1}\n".format(key, repr(value)))

  # aprun on Fucking Crays. mpirun everywhere else.
  file.write("aprun " if "NERSC_HOST" in environ else "mpirun ")
  file.write("-n {0} python {1} ".format(mppwidth, pyscript) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --{0}".format(key))
    else:             file.write(" --{0} {1}".format(key, value))
  file.write(" " + pickle + "\n")

def default_slurm( file, walltime = "05:45:00", mppwidth = 8, ppernode=8, queue = None,\
                   name = None, pyscript = None,  pickle = "job_pickle",\
                   outdir = None, **kwargs):
  """ Creates default slurm-script. Does not launch. 

      :Parameters:
        file
          File object to which to write. 
        walltime
          must be a string in "hh:mm:ss" format, or anything acceptable to the
          PBS implementation.
        mppwidth
          Number of processes (not processors) to use.
        ppernode
          Number of processes per node.
        queue
          Queue to use. If none will let slurm use default queue.
        name
          Name of the job.
        pyscript
          Python script to launch. Default: L{runme.py}
        pickle
          Fileame of job-tree pickle.
  """
  from os import environ
  from os.path import exists, join, abspath, dirname

  nnodes = mppwidth / ppernode + (0 if mppwidth % ppernode == 0 else 1)

  file.write("#! /bin/bash\n"\
             "#SBATCH --time={0}\n"\
             "#SBATCH --account=BES000\n"\
             "#SBATCH -N {1}\n"\
             "#SBATCH -n {2}\n".format(walltime, nnodes, mppwidth)) 
  if queue != None: file.write("#SBATCH -p {0}\n".format(queue))
  pbsdir = dirname(file.name)
  if name != None:
    file.write("#SBATCH -J {1} \n"\
               "#SBATCH -e \"{0}/err.{1}.%j\"\n"\
               "#SBATCH -o \"{0}/out.{1}.%j\"\n".format(pbsdir, name))
  else:
    file.write("#SBATCH -e \"{0}/err.%j\"\n"\
               "#SBATCH -o \"{0}/out.%j\"\n".format(pbsdir))
  if outdir != None: file.write("#SBATCH -D {0}\n".format(abspath(outdir)))

  file.write( "mpirun -np {0} numa_wrapper -ppn={1} python {2} "\
              .format(mppwidth, ppernode, pyscript) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --{0}".format(key))
    else:             file.write(" --{0} {1}".format(key, value))
  file.write(" " + pickle + "\n")

