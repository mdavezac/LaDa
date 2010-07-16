""" IPython functions and data. """
from . import jobs

class _Global(object):
  current = None
  """ Current dictionary. """

  pickle_filename = None
  """ Current pickle filename. """

  pickle_directory = None
  """ Directory of current pickle. """

def explore(self, arg):
  """ Starts exploration of a pickled jobdictionary. """
  from os.path import exists, split as splitpath, join, abspath
  from cPickle import load
  from .opt import open_exclusive
  current          = _Global().current
  pickle_filename  = _Global().pickle_filename
  pickle_directory = _Global().pickle_directory
  
  ip = self.api
  args = arg.split()
  if len(args) == 0:
    if current == None:
      print "No current jobs."
    elif pickle_filename == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Filename of jobdictionary: ", join(pickle_filename, pickle_directory)

  if isinstance(arg, jobs.JobDict):
    _Global.current = arg
    _Global.pickle_filename = None
    _Global.pickle_dictionary = None
    return

  def get_dict(filename, is_a_file = None):
    if is_a_file == None: 
      is_a_file = exists(filename)
      is_a_global = filename in ip.user_ns
      if is_a_global: is_a_global = isinstance(ip.user_ns[filename], jobs.JobDict)
      if is_a_global and is_a_file:
        print "Found both a file and a JobDict variable named %s.\n"\
              "Please use \"explore file %s\" or \"explore JobDict %s\"." \
              % (filename, filename, filename)
        return
      if not (is_a_file or is_a_global):
        print "Could not find either JobDict variable or a file with name %s." % (filename)
        return 
    else: is_a_global = not is_a_file

    if is_a_file:
      if not exists(filename): # checks the file exists.
        print "Could not find file %s." % (filename)
        return 
      with open_exclusive(filename, "r") as file: result = load(file)
      return result, abspath(splitpath(filename)[0]), splitpath(filename)[1]

    if is_a_global: # checks the variable exists.
      if not filename in ip.user_ns: 
        print "Could not find JobDict variable %s." % (filename)
        return
      elif not isinstance(ip.user_ns[filename], jobs.JobDict): # checks the file exists.
        print "Variable %s is not a JobDict object." % (filename)
        return
      return ip.user_ns[filename], None, None

  if len(args) == 1: 
    result = get_dict(args[0])
    if result != None:
      _Global.current = result[0]
      _Global.pickle_filename = result[1]
      _Global.pickle_directory = result[2]
    return 

def goto(self, arg):
  """ Moves current dictionary position and working directory (if appropriate). """
  from os import chdir
  from os.path import exists
  current          = _Global().current
  pickle_filename  = _Global().pickle_filename
  pickle_directory = _Global().pickle_directory
  if len(arg.split()) == 0:
    if current == None:
      print "No current jobs."
    elif pickle_filename == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Filename of jobdictionary: ", join(pickle_filename, pickle_directory)
  try: result = _Global().current[arg] 
  except KeyError: 
    print "Could not find %s in current directory." % (arg)
    return 

  if not (hasattr(result, "parent") and hasattr(result, "children")  and hasattr(result, "jobparams")):
    print "%s is a job parameter, not a directory."
    return
  _Global.current = result
  if pickle_directory != None: return
  chdir( join(_Global.pickle_directory, _Global.current.name[1:]) )
  return

def current_dictionary(self, arg): 
  """ Returns current dictionary. """
  return _Global.current
  

def _main():
  import IPython.ipapi
  ip = IPython.ipapi.get()
  ip.expose_magic("explore", explore)
  ip.expose_magic("current_dictionary", current_dictionary)
  ip.expose_magic("goto", goto)

_main()
