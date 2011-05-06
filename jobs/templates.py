""" Templates for job submissal. """
__docformat__ = "restructuredtext en"

def hopper_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = None, name = None, \
                pyscript = None, pickle = "job_pickle", outdir = None, header = None, **kwargs):
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
          Filename of job-tree pickle.
  """
  from os import environ
  from os.path import dirname, abspath

  pbsdir = abspath(dirname(file.name))
  file.write("#! /bin/bash\n#PBS -l walltime={0},mppwidth={1}\n".format(walltime, mppwidth))
  if name != None: 
    file.write("#PBS -N {1}\n"\
               "#PBS -e {0}/err.{1}.$PBS_JOBID\n"\
               "#PBS -o {0}/out.{1}.$PBS_JOBID\n".format(pbsdir, name))
  else:
    file.write("#PBS -e {0}/err.$PBS_JOBID\n"\
               "#PBS -o {0}/out.$PBS_JOBID\n".format(pbsdir))
  if queue != None: file.write("#PBS -q {0} \n".format(queue))
  file.write( "export XTPE_LINK_TYPE=dynamic\n"\
              "export CRAY_ROOTFS=DSL\n\n"\
              "module unload xt-libsci\n"\
              "module load python\n"\
              "module unload PrgEnv-pgi >& /dev/null\n"\
              "module load PrgEnv-gnu\n"\
              "module unload xt-shmem >& /dev/null\n"\
              "module load xt-mpich2\n\n" )

  if header != None: file.write("{0}\n".format(header))

  file.write("\nsource {0}_mpi/bin/preactivate"\
             "\nsource {0}_mpi/bin/activate"\
             "\nsource {0}_mpi/bin/postactivate"\
             "\nsource {0}/../preactivate"\
             "\nsource {0}/../postactivate\n\n".format(environ['VIRTUAL_ENV']))

  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd {0}\n".format(outdir))

  # aprun on fucking Cray prima donas. mpirun everywhere else.
  file.write("aprun" if environ.get('NERSC_HOST', 'other') == 'hopper2' else 'mpirun')
  file.write(" -n {0} python {1} ".format(mppwidth, pyscript) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --{0}".format(key))
    else:             file.write(" --{0} {1}".format(key, value))
  file.write(" " + pickle + "\n")

def default_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = None, name = None, \
                 pyscript = None, pickle = "job_pickle", outdir = None, header = None, **kwargs):
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
          Filename of job-tree pickle.
  """
  from os import environ
  from os.path import dirname, abspath

  pbsdir = abspath(dirname(file.name))
  file.write("#! /bin/bash\n#PBS -l walltime={0},mppwidth={1}\n".format(walltime, mppwidth))
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
  if header != None: file.write("{0}\n".format(header))
  if outdir == None: file.write("cd $PBS_O_WORKDIR\n")
  else: file.write("cd {0}\n".format(outdir))

  # aprun on fucking Cray prima donas. mpirun everywhere else.
  file.write("aprun" if environ.get('NERSC_HOST', 'other') == 'hopper2' else 'mpirun')
  file.write(" -n {0} python {1} ".format(mppwidth, pyscript) )
  for key, value in kwargs.items(): 
    if value == None: file.write(" --{0}".format(key))
    else:             file.write(" --{0} {1}".format(key, value))
  file.write(" " + pickle + "\n")

def default_slurm( file, walltime = "05:45:00", mppwidth = 8, ppernode=8, account = None,\
                   name = None, pyscript = None,  pickle = "job_pickle",\
                   outdir = None, partition = None, **kwargs):
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
        account
          Account to use. Defaults to BES000.
        name
          Name of the job.
        pyscript
          Python script to launch. Default: L{runme.py}
        pickle
          Filename of job-tree pickle.
        partition 
          Possible partition within which to launch jobs.
  """
  from os.path import abspath, dirname

  nnodes = mppwidth / ppernode + (0 if mppwidth % ppernode == 0 else 1)

  if account == None: account = "BES000"
  file.write("#! /bin/bash\n"\
             "#SBATCH --account={3}\n"\
             "#SBATCH --time={0}\n"\
             "#SBATCH -N {1}\n"\
             "#SBATCH -n {2}\n".format(walltime, nnodes, mppwidth, account)) 
  if partition != None: file.write("#SBATCH -p {0}\n".format(partition))
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
