debug_queue = None
""" No debug queue on cx1. """

qsub_exe = "qsub"
""" Qsub executable. """
          
default_pbs = { 'walltime': "00:55:00", 'nnodes': 1, 'ppn': 4}
""" Defaults parameters filling the pbs script. """

def pbs_string(**kwargs):
  """ Returns pbs script. """
  if 'name' in kwargs:
    kwargs['name'] = kwargs['name'][:min(len(kwargs['name']), 15)]
  return "#! /bin/bash\n"\
             "#PBS -e \"{err}\"\n"\
             "#PBS -o \"{out}\"\n"\
             "#PBS -m n\n"\
             "#PBS -N {name}\n"\
             "#PBS -l select={nnodes}:ncpus={n}\n"\
             "#PBS -l walltime={walltime}\n"\
             "#PBS -V \n\n"\
             "cd {directory}\n"\
             "python {scriptcommand}\n".format(**kwargs)

default_comm = { 'n': 2, 'ppn': default_pbs['ppn'] }
""" Default mpirun parameters. """

mpirun_exe = "mpirun {placement} -n {n} {program}"
""" Command-line to launch external mpi programs. """

def vasp_program(vasp):
  """ Signifies the vasp executable. 
  
      It is expected that two vasp executable exist, a *normal* vasp, and a one
      compiled for non-collinear calculations.

  """
  return "vasp-4.6-nc" if getattr(vasp, 'lsorbit', False) == True              \
         else "vasp-4.6"

def ipython_qstat(self, arg):
  """ Prints jobs of current user. """
  from subprocess import Popen, PIPE
  from IPython.utils.text import SList
  # get user jobs ids
  jobs   = SList(Popen(['qstat', '-f'], stdout=PIPE)                           \
                .communicate()[0].split('\n'))
  names  = [ u[u.find('=')+1:].lstrip().rstrip()                               \
             for u in jobs.grep('Job_Name') ]
  mpps   = [int(u[u.find('=')+1:]) for u in jobs.grep('Resource_List.ncpus')]
  states = [ u[u.find('=')+1:].lstrip().rstrip()                               \
             for u in jobs.grep('job_state') ]
  ids    = [u[u.find(':')+1:].lstrip().rstrip() for u in jobs.grep('Job Id')]
  return SList([ "{0:>10} {1:>4} {2:>3} -- {3}".format(id, mpp, state, name)   \
                 for id, mpp, state, name in zip(ids, mpps, states, names)])
