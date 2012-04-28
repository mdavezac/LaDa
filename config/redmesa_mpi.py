debug_queue = "queue", "inter"
""" How to select the debug queue. 

    First part of the tuple is the keyword argument to modify when calling
    the pbs job, and the second is its value.
"""

accounts = ["BES000"]
""" List of slurm or pbs accounts allowed for use. 

    This is used by ipython's %launch magic function. 
    It is not required for slurm systems. 
    If empty, then %launch will not have a queue option.
"""

qsub_exe = "sbatch"
""" Qsub executable. """
          
default_pbs = { 'account': accounts[0], 'walltime': "06:00:00", 'nnodes': 1, 'ppn': 8}
""" Defaults parameters filling the pbs script. """

pbs_string =  "#! /bin/bash\n"\
              "#SBATCH --account={account}\n"\
              "#SBATCH --time={walltime}\n"\
              "#SBATCH -N {nnodes}\n"\
              "#SBATCH -e \"{err}.%j\"\n"\
              "#SBATCH -o \"{out}.%j\"\n"\
              "#SBATCH -J {name}\n"\
              "#SBATCH -D {directory}\n\n"\
              "python {scriptcommand}\n"
""" Default slurm script. """

default_comm = { 'n': 2, 'ppn': default_pbs['ppn'] }
""" Default mpirun parameters. """

mpirun_exe = "mpirun -np {n} {placement} numa_wrapper -ppn={ppn} {program}"
""" Command-line to launch external mpi programs. """

def ipython_qstat(self, arg):
  """ squeue --user=`whoami` -o "%7i %.3C %3t  --   %50j" """
  from subprocess import Popen, PIPE
  from IPython.utils.text import SList
  from getpass import getuser

  # finds user name.
  whoami = getuser()
  squeue = Popen(["squeue", "--user=" + whoami, "-o", "\"%7i %.3C %3t    %j\""], stdout=PIPE)
  result = squeue.stdout.read().rstrip().split('\n')
  result = SList([u[1:-1] for u in result[1:]])
  return result.grep(str(arg[1:-1]))
