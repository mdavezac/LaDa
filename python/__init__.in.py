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

version_info = (@LaDa_VERSION_MAJOR@, @LaDa_VERSION_MINOR@)
""" Tuple containing version info. """
version = "{0[0]}.{0[1]}".format(version_info)
""" String containing version info. """
lada_with_mpi = @do_use_mpi@
""" If True, should load MPI stuff.

    If False, should try and avoid loading MPI related stuff.
"""

# reads stuff from global configuration files.
global_dict = {"ladamodules": __all__}
if "jobs" in __all__:
  from lada.jobs.templates import default_pbs, default_slurm
  from opt import cpus_per_node
  global_dict["default_pbs"]   = default_pbs
  global_dict["default_slurm"] = default_slurm
  global_dict["cpus_per_node"] = cpus_per_node
  del default_pbs
  del default_slurm
  del cpus_per_node

# first configuration files installed with lada.
from os.path import exists, expanduser, expandvars, dirname, join
from glob import iglob
from os import environ
for filename in iglob(join(join(dirname(__file__), "config"), "*.py")):
  local_dict = {}
  execfile(filename, global_dict, local_dict)
  locals().update((k, v) for k, v in local_dict.iteritems() if k[0] != '_')

# then configuration files installed in a global config directory.
if "LADA_CONFIG_DIR" in environ: 
  for filename in iglob(join(environ["LADA_CONFIG_DIR"], "*.py")):
    local_dict = {}
    execfile(filename, global_dict, local_dict)
    locals().update((k, v) for k, v in local_dict.iteritems() if k[0] != '_')

# then user configuration file.
if exists(expandvars(expanduser('~/.lada'))):
  local_dict = {}
  execfile(expandvars(expanduser('~/.lada')), global_dict, local_dict)
  locals().update((k, v) for k, v in local_dict.iteritems() if k[0] != '_')

# clean up namespace
del global_dict
del exists
del expanduser
del expandvars
del dirname
del join
del iglob
del environ
