""" Root of all lada python packages and modules. 

    Configuration variables exist here. However, they are added within separate
    files. Which files will depend upon the user.

       - Files located in the config sub-directory where lada is installed
       - Files located in one of the directories specified by :envvar:`LADA_CONFIG_DIR`
       - In the user configuration file ~/.lada

    The files are read in that order. Within a given directory, files are read
    alphabetically. Later files can override previous files. Finally, all and
    any variable which are declared within these files will end up at the root
    of :py:mod:`lada`. Be nice, don't pollute yourself.

    .. envvar:: LADA_CONFIG_DIR

       Environment variable specifying the path(s) to the configuration directories.
    
    For which variables can be set (and, if relevant, which file) see lada.config.
"""
__docformat__ = "restructuredtext en"
__all__ = ["error", @which_packages@]
import error
__all__ = [@which_packages@]

version_info = (@LaDa_VERSION_MAJOR@, @LaDa_VERSION_MINOR@)
""" Tuple containing version info. """
version = "{0[0]}.{0[1]}".format(version_info)
""" String containing version info. """


# reads stuff from global configuration files.
# doing it in a function makes it easier to keep the lada namespace clean.
def _config_files():
  from os.path import exists, expanduser, expandvars, dirname, join
  from glob import iglob
  from os import environ

  # dictionary with stuff we want defined when reading config files.
  global_dict = {"ladamodules": __all__}
  if "jobs" in __all__:
    from lada.jobs.templates import default_pbs, default_slurm
    global_dict["default_pbs"]   = default_pbs
    global_dict["default_slurm"] = default_slurm

  local_dict = {}
  # first configuration files installed with lada.
  for filename in iglob(join(join(dirname(__file__), "config"), "*.py")):
    execfile(filename, global_dict, local_dict)

  # then configuration files installed in a global config directory.
  if "LADA_CONFIG_DIR" in environ: 
    for filename in iglob(join(environ["LADA_CONFIG_DIR"], "*.py")):
      execfile(filename, global_dict, local_dict)

  # then user configuration file.
  if exists(expandvars(expanduser('~/.lada'))):
    execfile(expandvars(expanduser('~/.lada')), global_dict, local_dict)
  return local_dict

# does actual config call.
locals().update((k, v) for k, v in _config_files().iteritems() if k[0] != '_')

# clean up namespace
del _config_files
