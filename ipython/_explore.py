""" IPython explore function and completer. """

def _glob_job_pickles(ip, arg):
  """ Returns list of candidate files and directory. """
  from ..jobs import JobDict

  # nothing there yet.
  if len(arg) == 0:
    result = [ u + "/" for u in ip.magic("mglob dir:*") ]
    result.extend([ u for u in ip.magic("mglob \"cont:JobDict pickle\" *") ])
  # Contains directory.
  elif arg.find('/') != -1:
    orig_dir = getcwd()
    new_dir = arg[:-arg[::-1].find('/')-1]
    chdir (new_dir)
    result = [ join(new_dir, u) + "/" for u in ip.magic("mglob dir:*") ]
    result.extend([ join(new_dir, u) for u in ip.magic("mglob \"cont:JobDict pickle\" *") ])
    chdir(orig_dir)
  else: # not a directory.
    result = [u + "/" for u in ip.magic("mglob dir:%s*" % (arg))]
    result.extend([u for u in ip.magic("mglob \"cont:JobDict pickle\" %s*" % (arg))])
  return result


def _glob_ipy_user_ns(ip, arg):
  """ Returns list of candidate python variables. """
  from ..jobs import JobDict
  
  pyargs = arg.split('.')  
  if len(pyargs) == 0: return []
  dictionary, name = ip.user_ns, None
  for u in pyargs[:-1]:
    if u not in dictionary: return results
    if hasattr(dictionary[u], "__dict__"):
      dictionary = dictionary[u].__dict__
      if name == None: name = u
      else: name += "." + u
  names = [u for u in dictionary if u[0] != '_']
  if name == None: 
    result = [u + "." for u in names if not isinstance(dictionary[u], JobDict)]
    result.extend([ u for u in names if isinstance(dictionary[u], JobDict)])
  else: 
    result = [name + "." + u + "." for u in names if not isinstance(dictionary[u], JobDict)]
    result.extend([ name + "." + u for u in names if isinstance(dictionary[u], JobDict)])
    return result
  return result

def explore_completer(self, event):
  """ Completer for explore. """ 
  from os import getcwd, chdir
  from os.path import isdir, join
  import IPython
  from ..jobs import JobDict
  from . import _get_current_job_params

  ip = self.api
  current, path = _get_current_job_params(self, 0)

  args = event.symbol.split()
  if len(args) == 0: 
    result = _glob_job_pickles(ip, "")
    result.extend(["errors", "results", "file", "JobDict"])
    return result
  elif len(args) == 1: 
    result = _glob_job_pickles(ip, args[0])
    if args[0].find('/') == -1:
      result.extend(["errors", "results", "file", "JobDict"])
      result.extend(_glob_ipy_user_ns(ip, args[0]))
    return result
  raise IPython.ipapi.TryNext


def explore(self, arg):
  """ Starts exploration of a pickled job dictionary. 
  
      Usage: 
      The most standard form is to simply load a job dictionary. All other
      job-dictionary magic functions will then use it.
      
      >>> explore path/to/jobdictionary_pickle

      If you have created a jobdictionary directly (rather than save it to
      disk), you can also load it as

      >>> explore jobdict_variable 

      In case of conflict between a pathname and a variable name, you can use
      the more explicit version.

      >>> explore file jobdict
      >>> explore JobDict jobdict

      You can load a dictionary and filter out successfull or unsuccessfull runs. 
      To explore errors only, use:
     
      >>> explore errors in path/to/job_pickle

      To explore only successful results, use:

      >>> explore results in path/to/job_pickle
  """
  from os.path import exists, split as splitpath, abspath, isfile
  from cPickle import load
  from ..opt import open_exclusive
  from ..jobs import JobDict
  from . import _get_current_job_params
  ip = self.api
  current, path = _get_current_job_params(self, 0)
  ip.user_ns.pop("_lada_error", None)
  
  args = [u for u in arg.split() if u not in ["in", "with"]]
  if len(args) == 0:
    if current == None:
      print "No current jobs."
    elif path == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Path to job dictionary: ", path
    return

  if isinstance(arg, JobDict):
    ip.user_ns["current_jobdict"] = arg
    ip.user_ns["current_jobdict_path"] = None
    return

  def get_dict(filename, is_a_file = None):
    if is_a_file == None: 
      is_a_file = exists(filename)
      var, dvar = None, ip.user_ns
      for name in  filename.split('.'):
        try: var, dvar = dvar[name], dvar[name].__dict__ 
        except: var = None; break
      is_a_variable = var != None
      if is_a_variable and is_a_file:
        raise RuntimeError(\
                "Found both a file and a JobDict variable named %s.\n"\
                "Please use \"explore file %s\" or \"explore JobDict %s\"." \
                % (filename, filename, filename) )
      if not (is_a_file or is_a_variable):
        raise RuntimeError(\
            "Could not find either JobDict variable or a file with name %s." % (filename) )
    else: is_a_variable = not is_a_file

    if is_a_file:
      if not exists(filename): # checks the file exists.
        raise RuntimeError("Could not find file %s." % (filename))
      if not isfile(filename):
        raise RuntimeError("%s is not a file." % (filename))
      with open_exclusive(filename, "r") as file: result = load(file)
      return result, abspath(filename)

    if is_a_variable: # checks the variable exists.
      var, dvar = None, ip.user_ns
      for name in  filename.split('.'):
        try: var, dvar = dvar[name], dvar[name].__dict__ 
        except:
          raise RuntimeError("Could not find variable name %s." % (filename))
      if not isinstance(var, JobDict): # checks the file exists.
        raise RuntimeError("Variable %s is not a JobDict object." % (filename))
      return var, None

  def _impl():
    """ Implementation of explore. """
    from os.path import join, exists
    from copy import deepcopy
    # filters some results, depending on success or failure.
    if args[0] == "errors" or args[0] == "results":
      # other arguments should be dictionary.
      string = ""
      for i in args[1:]: string += " " + i
      explore(self, string)
      current, path = _get_current_job_params(self, 0)
      if current == None:
        raise RuntimeError("Error encountered while trying to explore dictionary." )
      if path == None:
        raise RuntimeError("No path set for current job dictionary.\n"\
                           "Please set current_dictionary_path to correct value.\n")
      ip.user_ns["current_jobdict"] = deepcopy(ip.user_ns["current_jobdict"])
      rootdir = splitpath(path)[0]
      # now marks (un)successful runs.
      which = (lambda x: not x) if args[0] == "results" else (lambda x: x)
      for job, d in ip.user_ns["current_jobdict"].walk_through():
        if job.is_tagged: continue
        directory = join(rootdir, d)
        if which(job.functional.Extract(directory).success): job.tag()
        else: job.untag()
    elif args[0] == "jobs": 
      string = ""
      for i in args[1:]: string += " " + i
      explore(self, string)
      current, path = _get_current_job_params(self, 0)
    elif len(args) == 1:                          current, path = get_dict(args[0])
    elif len(args) == 2 and args[0] == "file":    current, path = get_dict(args[0], True)
    elif len(args) == 2 and args[0] == "JobDict": current, path = get_dict(args[0], False)
    else: raise RuntimeError("Calling explore with arguments %s is invalid." % (arg))
    if current != None: 
      ip.user_ns["current_jobdict"] = current
      ip.user_ns["current_jobdict_path"] = path

  try: _impl()
  except RuntimeError as e:
    ip.user_ns["_lada_error"] = e
    print e
  else:
    # remove any prior iterator stuff.
    ip.user_ns.pop("_lada_subjob_iterator", None)
    ip.user_ns.pop("_lada_subjob_iterated", None)

