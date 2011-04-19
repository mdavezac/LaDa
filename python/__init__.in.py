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
lada_with_slurm = 'SNLCLUSTER' in environ 
""" If True use slurm as ressource manager, else use openpbs. """
queues = []
""" List of slurm or pbs queues allowed for use. 

    This is used by ipython's %launch magic function. 
    It is not required for slurm systems. 
    If empty, then %launch will not have a queue option.
"""

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

# reads stuff from input file
from os.path import exists, expanduser, expandvars
path = expandvars(expanduser('~/.lada'))
if exists(path):
  from opt import read_input
  with open(path, 'r') as file: string = file.read()
  global_dict, local_dict = {}, {}
  exec(string, global_dict, local_dict)
  locals().update(local_dict)
