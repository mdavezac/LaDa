""" Lamarck-Darwin extension library.
    =================================


    Getting the source, installing, and compiling LaDa:
    ---------------------------------------------------
    Getting there.

    IPython interface for high-thoughput calculations:
    --------------------------------------------------

    This is an interface to monitor highthoughput calculations. It allows for
    launching job-dictionaries, rapid access to jobs which did not complete
    successfully, as well as output collection (e.g. total-energy,
    eigenvalues,...) across all calculation in a job-dictionary.  

    Please see `lada.ipython` for a more in-depth description.


    Creating a job-dictionary:
    --------------------------

    In the works.
    Please see `lada.jobs` for a more in-depth description.


    Constructing and manipulating crystal-structures:
    -------------------------------------------------

    In the works.
    Please see `lada.crystal` for a more in-depth description.

    Interacting with VASP:
    ----------------------

    In the works.
    Please see `lada.vasp` for a more in-depth description.

    Setting defaults:
    -----------------

    All global variables below can be changed within "$HOME/.lada". This file
    is read only once, when first importing lada into python.  
"""
__docformat__ = "restructuredtext en"
__all__ = [@which_packages@]
from os import environ

version_info = (@LaDa_VERSION_MAJOR@, @LaDa_VERSION_MINOR@)
""" Tuple containing version info. """
version = "{0[0]}.{0[1]}".format(version_info)
""" String containing version info. """
lada_with_mpi = @do_use_mpi@
""" If True, should load MPI stuff.

    If False, should try and avoid loading MPI related stuff.
"""
# right now just check for redmesa or redrock. 

genpot_library = "libgenpot.so"
""" Default genpot library. 

    The value for the default can be overriden by ~/.lada in the code below.
"""

escan_library = "libpescan.so"
""" Default escan library. 

    The value for the default can be overriden by ~/.lada in the code below.
"""

vasp_library = "libvasp.so"
""" Default vasp library. 

    The value for the default can be overriden by ~/.lada in the code below.
"""


# The variables defined below are only needed for the ipython interface.
if @with_ipython@:
  from os import environ
  from .jobs.templates import default_pbs, default_slurm

  default_walltime = "06:00:00"
  """ Default walltime when launching jobs. """

  lada_with_slurm = 'SNLCLUSTER' in environ 
  """ If True use slurm as ressource manager, else use openpbs. """
  queues = []
  """ List of slurm or pbs queues allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """
  accounts = []
  """ List of slurm or pbs accounts allowed for use. 

      This is used by ipython's %launch magic function. 
      It is not required for slurm systems. 
      If empty, then %launch will not have a queue option.
  """

  template_pbs = default_pbs
  """ Template pbs script to use. Depends on machine. """
  mpirun_exe = "mpirun"
  """ Name of mpirun-like executable. """

  debug_queue = "queue", "debug"
  """ How to select the debug queue. 

      First part of the tuple is the keyword argument to modify when calling
      the pbs job, and the second is its value.
  """
  qsub_exe = "qsub"
  """ Qsub executable. """
  resource_string = "nodes={1}:ppn={2}" 
  """ Format string to specify computational resources. 
      
      The first argument is total number of processes, the second the number of
      nodes itself, the third the number of processes per node.
  """

  if environ.get("SNLCLUSTER", "none") in ["redrock", "redmesa"]: 
    template_pbs = default_slurm
    debug_queue = "queue", "inter"
    accounts = ["BES000"]
    qsub_exe = "sbatch"
    resource_string = "-N {1}"
  elif environ.get("NERSC_HOST", "none") == "hopper":
    mpirun_exe = "aprun"
    queues = "debug", "regular", "low", "premimum"
    resource_string = "mppwidth={0}"
  elif environ.get("NERSC_HOST", "none") == "carver":
    queues = "debug", "regular", "low"
    resource_string = "nodes={1}:ppn=8"


# reads stuff from input file
from os.path import exists, expanduser, expandvars
if exists(expandvars(expanduser('~/.lada'))):
  from opt import read_input
  with open(expandvars(expanduser('~/.lada')), 'r') as file: string = file.read()
  global_dict, local_dict = {}, {}
  exec(string, global_dict, local_dict)
  locals().update(local_dict)
