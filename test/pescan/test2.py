from os.path import join, exists
from os import makedirs
from shutil import rmtree, copy
from boost.mpi import world
from lada.escan import call_escan as escan_as_library
from lada.opt.changedir import Changedir
from lada.opt import Redirect

# pseudo files
pseudo_files = [ "maskr", "vq.Ge", "vq.Si", "vq.SiGe.Ge", "vq.SiGe.Si", "vwr.pso" ]

# directories where jobs are at.
jobs = ["VBM", "Gamma",  "X",  "L",  "W1", "W2"]
testdir = "test_input"
workdir = "work"
if world.rank == 0 and (not exists(workdir)): makedirs(workdir)

# creates local comm split into N groups.
N = 2
color = world.rank % N
local_comm = world.split(color)

# launch pescan for different jobs.
for i, dir in enumerate(jobs):
  # splits job according to color.
  if i % N != color: continue

  # creates working directory with all input files.
  workhere = join(workdir, dir)
  if local_comm.rank == 0:
    if exists(workhere): rmtree(workhere) # deletes everything if already there.
    makedirs(workhere)
    # symlink potential files
    testhere = join(testdir, "pseudos")
    for file in pseudo_files: copy(join(testhere, file), workhere)
    # symlinks input files.
    testhere = join(testdir, dir)
    for file in ["atom_config", "pot_input", "escan_input"]: 
      copy(join(testhere, file), workhere)

  # sync all coms. 
  local_comm.barrier()
  # changes current directory to working directory.
  with Changedir(workhere) as current_dir:
    # Redirects fortran output and error.
    stdout = "stdout"
    stderr = "stderr"
    if local_comm.rank != 0: 
      stdout = "%s.%i" % (stdout, local_comm.rank)
      stderr = "%s.%i" % (stderr, local_comm.rank)
    with Redirect(Redirect.fortran.output, stdout) as stdout:
      with Redirect(Redirect.fortran.error, stderr) as stderr:
        # finally calls escan.
        escan_as_library( local_comm, \
                          atom="atom_config",\
                          pot="pot_input", \
                          escan="escan_input" )
  
