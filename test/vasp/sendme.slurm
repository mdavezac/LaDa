#! /bin/bash
#SBATCH --time=24:00:00             # WALLTIME
#SBATCH -N  6                       # Number of nodes
#SBATCH --job-name vaspfarm              # Name of job
#SBATCH -e err
#SBATCH -o out
nodes=$SLURM_JOB_NUM_NODES          # Number of nodes
cores=16                            # Number MPI processes to start on each node

if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi
# setup bash environment ...
module load libraries/intel-mkl-11.1.046 
if [ -n "$LD_LIBRARY_PATH" ] ; then
  LD_LIBRARY_PATH="/home/mdadcast/usr/lib:$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH
else
  LD_LIBRARY_PATH="/home/mdadcast/usr/lib"
  export LD_LIBRARY_PATH
fi
export PATH=/bin:$PATH
# ... so as to setup python environment
source /home/mdadcast/virtualenv/vasp/bin/activate

# executes python script.
echo $LD_LIBRARY_PATH
rm -rf work
mpiexec -np $(($nodes*$cores)) python test.py
