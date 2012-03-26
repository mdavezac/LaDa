mpirun_exe = "mpirun -n {n} {program}"
""" Command-line to launch external mpi programs. """

default_comm = {'n': 2, 'ppn': 4}
""" Default communication directory. 

    should contain all key-value pairs used in :py:data:`mpirun_exe`.
"""

# pbs/slurm related stuff.
queues = ()
""" List of slurm or pbs queues allowed for use. 

    This is used by ipython's %launch magic function. 
    It is not required for slurm systems. 
    If empty, then %launch will not have a queue option.
"""
accounts = ['BES000']
""" List of slurm or pbs accounts allowed for use. 

    This is used by ipython's %launch magic function. 
    It is not required for slurm systems. 
    If empty, then %launch will not have a queue option.
"""

debug_queue = "queue", "debug"
""" How to select the debug queue. 

    First part of the tuple is the keyword argument to modify when calling
    the pbs job, and the second is its value.
"""
qsub_exe = "sbatch"
""" Qsub/sbatch executable. """

default_pbs = {'account': accounts[0], 'walltime': "06:00:00", 'nnodes': 1 }
""" Defaults parameters filling the pbs script. """
pbs_string =  "#! /bin/bash/\n"\
              "#SBATCH --account={account}\n"\
              "#SBATCH --time={walltime}\n"\
              "#SBATCH -N={nnodes}\n"\
              "#SBATCH -e={err}\n"\
              "#SBATCH -o={out}\n"\
              "#SBATCH -J={name}\n"\
              "#SBATCH -D={directory}\n\n"\
              "python {scriptcommand}\n"
""" Default pbs/slurm script. """
