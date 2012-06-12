debug_queue = None
""" No debug queue on cx1. """

qsub_exe = "qsub"
""" Qsub executable. """
          
default_pbs = { 'walltime': "00:55:00", 'nnodes': 1, 'ppn': 4}
""" Defaults parameters filling the pbs script. """

pbs_string =  "#! /bin/bash\n"\
              "#PBS -e \"{err}.%j\"\n"\
              "#PBS -o \"{out}.%j\"\n"\
              "#PBS -m n\n"\
              "#PBS -N {name}\n"\
              "#PBS -l select={nnodes}:ncpus={nprocs}\n"\
              "#PBS -l walltime={walltime}\n"\
              "#PBS -V \n\n"\
              "cd {directory}\n"\
              "python {scriptcommand}\n"
""" Default slurm script. """

default_comm = { 'n': 2, 'ppn': default_pbs['ppn'] }
""" Default mpirun parameters. """

mpirun_exe = "mpirun {placement} -n {n} {program}"
""" Command-line to launch external mpi programs. """

# def placement(communicator=None):
#   """ Placement string for MPI processes. 
#   
#       Should return an empty string if the communicator is None.
#       This version works for openmpi.
#   """
#   print 'WTF'
#   if communicator is None: return ""
#   if len(getattr(communicator, 'machines', {})) == 0: return ""
#   return "-f {0}".format(communicator.nodefile())

def vasp_program(self):
  """ Signifies the vasp executable. 
  
      It is expected that two vasp executable exist, a *normal* vasp, and a one
      compiled for non-collinear calculations.

  """
  return "vasp-4.6-nc" if getattr(self, 'lsorbit', False) == True else "vasp-4.6"

def ipython_qstat(self, arg):
  """ Prints jobs of current user. """
  from subprocess import Popen, PIPE
  from IPython.utils.text import SList
  # get user jobs ids
  jobs   = SList(Popen(['qstat', '-f'], stdout=PIPE).communicate()[0].split('\n'))
  names  = [u[u.find('=')+1:].lstrip().rstrip() for u in jobs.grep('Job_Name')]
  mpps   = [int(u[u.find('=')+1:]) for u in jobs.grep('Resource_List.ncpus')]
  states = [u[u.find('=')+1:].lstrip().rstrip() for u in jobs.grep('job_state')]
  ids    = [u[u.find(':')+1:].lstrip().rstrip() for u in jobs.grep('Job Id')]
  return SList([ "{0:>10} {1:>4} {2:>3} -- {3}".format(id, mpp, state, name) \
                 for id, mpp, state, name in zip(ids, mpps, states, names)])
