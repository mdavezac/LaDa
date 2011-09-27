""" Templates for job submissal. """
__docformat__ = "restructuredtext en"

def default_pbs( file, walltime=None, mppwidth=8, ppernode=None, queue=None, \
                 account=None, name=None,  pyscript=None, pickle="job_pickle", \
                 outdir=None, external=False, memlim=None, **kwargs):
  """ Creates pbs-script for carver at nersc. Does not launch. 

      :Parameters:
        file
          File object to which to write. 
        walltime
         must be a string in "hh:mm:ss" format, or anything
         acceptable to the PBS implementation. If None, defaults to `lada.default_walltime`.
        mppwidth
          Number of processes (not processors) to use.
        ppernode 
          Number of processes per node. If None, uses `lada.cpus_per_node` in lada.
        queue
          Queue to use
        account
          Account under which to charge computer time, if any.
        name
          Name of the job.
        pyscript
          Python script to launch. Default: L{runme.py}
        pickle
          Filename of job-dictionary pickle.
        outdir 
          Directory where to start calculations.
        external
          Whether computations are launched as libraries or external jobs.
          If external, then the main script is launched without MPI.
        memlim
          Memory limits to impose with ulimit. If 0, don't impose limits.
          If negative, guess from /proc/meminfo and number of processes per node.
        kwargs
          Further arguments are appended at the end of mpirun.
  """
  from os.path import dirname, abspath
  from .. import mpirun_exe, resource_string, default_walltime, cpus_per_node

  if walltime == None: walltime = default_walltime 
  if ppernode == None: ppernode = cpus_per_node

  pbsdir = abspath(dirname(file.name))
  nnodes = mppwidth//8 if mppwidth % 8 == 0 else mppwidth//8 + 1
  file.write("#! /bin/bash\n#PBS -l walltime={0}\n".format(walltime))
  file.write("#PBS -l " + resource_string.format(mppwidth, nnodes, ppernode) + "\n")
  if name != None: 
    file.write("#PBS -N {1}\n\n"\
               "#PBS -V\n"\
               "#PBS -e {0}/err.{1}.$PBS_JOBID\n"\
               "#PBS -o {0}/out.{1}.$PBS_JOBID\n".format(pbsdir, name))
  else:
    file.write("#PBS -V\n"\
               "#PBS -e {0}/err.$PBS_JOBID\n"\
               "#PBS -o {0}/out.$PBS_JOBID\n".format(pbsdir))
  if queue != None: file.write("#PBS -q {0} \n".format(queue))
  if account != None: file.write("#PBS -A {0} \n".format(account))
  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd {0}\n".format(outdir))
  if memlim < 0:
    file.write( "ulimit -v `python -c \"from lada.opt import total_memory; total_memory() / {0}\n\"`"\
                .format(ppernode) )
  elif memlim > 0:
    file.write( "ulimit -v {0}".format(memlim) )

  # aprun on fucking Cray prima donas. mpirun everywhere else.
  if external:
    file.write("python {0} --nprocs {1}".format(pyscript, mppwidth))
  else:
    file.write(mpirun_exe.format(mppwidth, "python {0}".format(pyscript)) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --{0}".format(key))
    else:             file.write(" --{0} {1}".format(key, value))
  file.write(" " + pickle + "\n")

def default_slurm( file, walltime = "05:45:00", mppwidth = 8, ppernode=None, account=None,\
                   name=None, pyscript=None,  pickle="job_pickle", queue=None, \
                   outdir=None, external=False, memlim=None, **kwargs):
  """ Creates default slurm-script. Does not launch. 

      :Parameters:
        file
          File object to which to write. 
        walltime
          must be a string in "hh:mm:ss" format, or anything acceptable to the
          PBS implementation.
        mppwidth
          umber of processes (not processors) to use.
        ppernode
          Number of processes per node. If None, uses `lada.cpus_per_node`.
        account
          Account to use. Defaults to BES000.
        name
          Name of the job.
        pyscript
          Python script to launch. Default: L{runme.py}
        pickle
          Filename of job-tree pickle.
        queue 
          Possible queue/partition within which to launch jobs.
        outdir 
          Working directory where to start job.
        memlim
          Memory limits to impose with ulimit. If None, don't impose limits.
          If "guess", guess from /proc/meminfo and number of processes per node.
        kwargs
          Further arguments are appended at the end of mpirun.
  """
  from os.path import abspath, dirname
  from .. import mpirun_exe, resource_string, default_walltime, cpus_per_node

  if walltime == None: walltime = default_walltime 
  if ppernode == None: ppernode = cpus_per_node
  nnodes = mppwidth // ppernode + (0 if mppwidth % ppernode == 0 else 1)

  if account == None: account = "BES000"
  file.write("#! /bin/bash\n"\
             "#SBATCH --account={1}\n"\
             "#SBATCH --time={0}\n".format(walltime, account)) 
  file.write("#SBATCH " + resource_string.format(mppwidth, nnodes, ppernode) + "\n")
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
  if memlim == "guess":
    file.write( "ulimit -v `python -c \"from lada.opt import total_memory; total_memory() / {0}\n\"`"\
                .format(ppernode) )
  elif memlim != None: 
    file.write( "ulimit -v {0}".format(memlim) )

  if external:
    file.write("python {0} --nprocs {1}".format(pyscript, mppwidth))
  else:
    file.write(mpirun_exe.format(nprocs=mppwidth, ppernode=ppernode, program="python {0}".format(pyscript)) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --{0}".format(key))
    else:             file.write(" --{0} {1}".format(key, value))
  file.write(" " + pickle + "\n")
