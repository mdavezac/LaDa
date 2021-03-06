""" Root of all pylada python packages and modules. 

    Configuration variables exist here. However, they are added within separate
    files. Which files will depend upon the user.

       - Files located in the config sub-directory where pylada is installed
       - Files located in one of the directories specified by :envvar:`PYLADA_CONFIG_DIR`
       - In the user configuration file ~/.pylada

    The files are read in that order. Within a given directory, files are read
    alphabetically. Later files can override previous files. Finally, all and
    any variable which are declared within these files will end up at the root
    of :py:mod:`pylada`. Be nice, don't pollute yourself.

    .. envvar:: PYLADA_CONFIG_DIR

       Environment variable specifying the path(s) to the configuration directories.
    
    For which variables can be set (and, if relevant, which file) see pylada.config.
"""
__docformat__ = "restructuredtext en"
__all__ = ["error", @which_packages@]
import error
import physics 

version_info = (@Pylada_VERSION_MAJOR@, @Pylada_VERSION_MINOR@)
""" Tuple containing version info. """
version = "{0[0]}.{0[1]}".format(version_info)
""" String containing version info. """


# reads stuff from global configuration files.
# doing it in a function makes it easier to keep the pylada namespace clean.
def _config_files(dointeractive=False):
  from os.path import exists, expanduser, expandvars, dirname, join
  from glob import iglob
  from os import environ

  # pattern to distinguish files to run only in interactive mode.
  # these files are loaded by the pylada-ipython extension itself.
  pattern = "*.py" if not dointeractive else "ipy_*.py"
  # dictionary with stuff we want defined when reading config files.
  global_dict = {"pyladamodules": __all__}
  local_dict = {}
  # first configuration files installed with pylada.
  for filename in iglob(join(join(dirname(__file__), "config"), pattern)):
    if dointeractive == False and filename[:4] == 'ipy_': continue
    execfile(filename, global_dict, local_dict)

  # then configuration files installed in a global config directory.
  if "PYLADA_CONFIG_DIR" in environ: 
    for directory in environ["PYLADA_CONFIG_DIR"].split(':'):
      for filename in iglob(join(directory, pattern)):
        if dointeractive == False and filename[:4] == 'ipy_': continue
        execfile(filename, global_dict, local_dict)

  # then user configuration file.
  if exists(expandvars(expanduser('~/.pylada'))):
    execfile(expandvars(expanduser('~/.pylada')), global_dict, local_dict)
  return local_dict

# does actual config call.
locals().update((k, v) for k, v in _config_files().iteritems() if k[0] != '_')
