""" Templates for job submissal. """

def default_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = "regular", name = None, \
                 pyvirt = "vasp.4.6.32", pbspools = 1, npbs = 0, procpools = 1, \
                 pyscript = None, pickle = "job_pickle", **kwargs ):
  from os.path import exists

  if name == None:
    name = "automatic" if pbspools <= 1 else "automatic_%i" % (npbs)

  if pyscript == None: # just copy standard script.
    pyscript = __file__.replace("templates.py", "runme.py")
    if pyscript[-3:] == "pyc": pyscript = pyscript[:-1]
  if isinstance(pyscript, str) and exists(pyscript): # just copy standard script.
    with open(pyscript, "r") as runme_file: pyscript = runme_file.readlines()

  file.write("#! /bin/bash\n")
  file.write("#PBS -l walltime=%s,mppwidth=%i\n" % (walltime, mppwidth) )
  file.write("#PBS -q " + queue + "\n")
  file.write("#PBS -e err\n")
  file.write("#PBS -o out\n")
  file.write("#PBS -N %s \n\n" % (name))

  file.write("cd $PBS_O_WORKDIR\n")
  file.write("workon %s \n" % (pyvirt) )

  file.write("\n\ncat << EOF >> runme.py\n")
  if isinstance(pyscript, str): file.write(pyscript)
  else:                         file.writelines(pyscript)
  file.write("\nEOF\n\n")

  file.write("mpirun -n %i python runme.py " % (mppwidth) )
  if pbspools > 1: file.write("--pbspools %i --npbs %i" % (pbspools, npbs) )
  if procpools > 1: file.write("--procpools %i" % (procpools) )
  file.write(" " + pickle + "\n")

def hopper_pbs( file, walltime = "05:45:00", mppwidth = 8, queue = "regular", name = None, \
                pyvirt = "vasp.4.6.32", pbspools = 1, npbs = 0, procpools = 1, \
                pyscript = None, pickle = "job_pickle", **kwargs ):

  if name == None:
    name = "automatic" if pbspools <= 1 else "automatic_%i" % (npbs)

  if pyscript == None: # just copy standard script.
    pyscript = __file__.replace("templates.py", "runme.py")
    if pyscript[-3:] == "pyc": pyscript = pyscript[:-1]
  if isinstance(pyscript, str) and exists(pyscript): # just copy standard script.
    with open(pyscript, "r") as runme_file: pyscript = runme_file.readlines()

  file.write("#! /bin/bash\n")
  file.write("#PBS -l walltime=%s,mppwidth=%i\n" % (walltime, mppwidth) )
  file.write("#PBS -q " + queue + "\n")
  file.write("#PBS -e err\n")
  file.write("#PBS -o out\n")
  file.write("#PBS -N %s \n\n" % (name))
  file.write("cd $PBS_O_WORKDIR\n")
  file.write("export CRAY_ROOTFS=DSL\n")
  file.write("workon %s \n" % (pyvirt) )
  file.write("module swap PrgEnv-pgi PrgEnv-gnu\n")
  file.write("module unload xt-libsci\n")
  file.write("module load fftw/2.1.5.2 pgi\n\n")

  file.write("\n\ncat << EOF >> runme.py\n")
  file.writelines(pyscript)
  file.write("\nEOF\n\n")
  file.write("aprun -n %i python runme.py " % (mppwidth) )
  if pbspools > 1: file.write("--pbspools %i --npbs %i" % (pbspools, npbs) )
  if procpools > 1: file.write("--procpools %i" % (procpools) )
  file.write(" " + pickle + "\n")

