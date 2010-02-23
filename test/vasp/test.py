from os.path import join, exists
from os import makedirs, getcwd
from shutil import rmtree, copy
from boost.mpi import world
from lada.vasp import call_vasp as vasp_as_library, version as vasp_version
from lada.opt.changedir import Changedir
from lada.opt import Redirect

request_vasp = "5.2.2"
assert vasp_version == request_vasp, \
       "Wrong vasp version: is %s, you requested %s." % (vasp_version, request_vasp)

# Creates working directory
workdir = "work"
assert not exists(workdir), "Please delete %s " % (workdir)
world.barrier()
if world.rank == 0:
  makedirs(workdir)

# pseudo files
input_dir = join(getcwd(), "input_files")
pseudo_files = [ "KPOINTS", "POSCAR", "POTCAR" ]

# Jobs
scale = 3e0
jobs = [ (Ueff / scale, Ueff*(1e0-0.5/scale), head) \
         for Ueff in [1e0, 0.5, 5] \
         for head in range(2,5)] 


# creates local comm split into N groups.
N = 6
color = world.rank % N
local_comm = world.split(color)

# launch pescan for different jobs.
for i, (Ueff, Vnlep, head) in enumerate(jobs):
  # splits job according to color.
  if i % N != color: continue

  # working directory with all input files.
  workhere = join(workdir, "U:%f-V:%f-head:%i" % (Ueff, Vnlep, head))

  # *********************************     only first comm writes
  if local_comm.rank == 0: 
    if exists(workhere): rmtree(workhere) # deletes everything if already there.
    makedirs(workhere)
    # symlink potential files
    for file in pseudo_files: copy(join(input_dir, file), workhere)
    # create INCAR
    with open(join(workhere, "INCAR"), "w") as incar:
      # writes head.
      with open(join(input_dir, "head%i" %(head)), "r") as headfile:
        incar.write(headfile.read()) 
      # writes the rest of it.
      nlep_string = \
          """# NLEP parameters                                              \n"""\
          """LDAU     = .TRUE.                                              \n"""\
          """LDAUPRINT= 1                                                   \n"""\
          """LDAUTYPE = 1                                                   \n"""\
          """                                                               \n"""\
          """LDUL1    =  -1  2        1                                     \n"""\
          """LDUU1    = 0.0  %f 1.255  ! Co d and O p orbital               \n"""\
          """LDUJ1    = 0.0  0.0      0.0                                   \n"""\
          """LDUO1    = 1    2        2                                     \n"""\
          """                                                               \n"""\
          """LDUL2    = -1    2      1                                      \n"""\
          """LDUU2    =  0.0 %f 6.493   ! Co d  nad O p orbital             \n"""\
          """LDUJ2    =  0.0 0.0     0.0                                    \n"""\
          """LDUO2    =  1   1       1         ! 1 means LDAU, 2 means NLEP \n"""\
          """                                                               \n"""\
          """LDUL3    =  2    2      -1                                     \n"""\
          """LDUU3    =  6.0 3.0     0.0   ! d orbital                      \n"""\
          """LDUJ3    =  0.0 1.0     0.0                                    \n"""\
          """LDUO3    =  1  1  1         ! 1 means LDAU, 2 means NLEP       \n"""\
          %  (Vnlep, Ueff)
      incar.write(nlep_string)
  # *********************************     only first comm writes

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
        vasp_as_library( local_comm )
  
