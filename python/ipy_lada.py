""" IPython functions and data. """
from . import jobs

def explore(self, arg):
  """ Starts exploration of a pickled jobdictionary. """
  from os.path import exists, split as splitpath, join, abspath
  from cPickle import load
  from .opt import open_exclusive
  ip = self.api
  current          = ip.user_ns.pop("current_jobdict", None)
  pickle_filename  = ip.user_ns.pop("current_jobdict_filename", None)
  pickle_directory = ip.user_ns.pop("current_jobdict_directory", None)
  
  args = arg.split()
  if len(args) == 0:
    if current == None:
      print "No current jobs."
    elif pickle_filename == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Path to job dictionary: ", join(pickle_filename, pickle_directory)

  if isinstance(arg, jobs.JobDict):
    ip.user_ns["current_jobdict"] = arg
    ip.user_ns["current_jobdict_filename"] = None
    ip.user_ns["current_jobdict_directory"] = None
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
      ip.user_ns["current_jobdict"] = result[0]
      ip.user_ns["current_jobdict_directory"] = result[1]
      ip.user_ns["current_jobdict_filename"] = result[2]
    return 
  elif len(args) == 2 and args[0] == "file":
    result = get_dict(args[0], True)
    if result != None:
      ip.user_ns["current_jobdict"] = result[0]
      ip.user_ns["current_jobdict_directory"] = result[1]
      ip.user_ns["current_jobdict_filename"] = result[2]
    return 
  elif len(args) == 2 and args[0] == "JobDict":
    result = get_dict(args[0], False)
    if result != None:
      ip.user_ns["current_jobdict"] = result[0]
      ip.user_ns["current_jobdict_directory"] = result[1]
      ip.user_ns["current_jobdict_filename"] = result[2]
    return 
  else: 
    print "Calling explore with arguments %s is invalid." % (arg)
    return

def goto(self, arg):
  """ Moves current dictionary position and working directory (if appropriate). """
  from os import chdir
  from os.path import exists, join
  ip = self.api
  current = ip.user_ns["current_jobdict"] if "current_jobdict" in ip.user_ns else None
  pickle_filename = ip.user_ns["current_jobdict_filename"] \
                    if "current_jobdict_filename" in ip.user_ns else None
  pickle_directory = ip.user_ns["current_jobdict_directory"] \
                     if "current_jobdict_directory" in ip.user_ns else None
    
  if current == None: 
    print "No current jobs."
    return
  if len(arg.split()) == 0:
    if pickle_filename == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Filename of jobdictionary: ", join(pickle_filename, pickle_directory)
    return
  arg = arg.split()[0]
  try: result = current[arg] 
  except KeyError as e: 
    print e
    return 

  if not (hasattr(result, "parent") and hasattr(result, "children")  and hasattr(result, "jobparams")):
    print "%s is a job parameter, not a directory."
    return
  ip.user_ns["current_jobdict"] = result
  current = ip.user_ns["current_jobdict"]
  if pickle_directory == None: return
  dir = join(pickle_directory, current.name[1:]) 
  if exists(dir): chdir(dir)
  return

def listjobs(self, arg):
  """ Lists subjobs. """
  ip = self.api
  if "current_jobdict" not in ip.user_ns: 
    print "No current job defined. Please first load a job with \"explore\"."
    return
  current = ip.user_ns["current_jobdict"] 
  if len(arg) != 0:
    try: subdict = ip.user_ns["current_jobdict"][arg] 
    except KeyError:
      print "%s is not a valid jobname of current job dictionary." % (arg)
      return
    current = current[arg]
    if not hasattr(current, "children"):  
      print "%s is not a valid jobname of current job dictionary." % (arg)
      return
  if len(current.children) == 0: return
  string = ""
  lines = ""
  for j in current.children.keys():
    if len(string+j) > 60:
      if len(lines) != 0: lines += "\n" + string
      else: lines = string
      string = ""
    string += j + " "

  if len(lines) == 0: print string
  else: print lines + "\n" + string
 
def goto_completer(self, event):
  import IPython
  ip = self.api
  if "current_jobdict" not in ip.user_ns: raise IPython.ipapi.TryNext
  elif '/' in event.symbol:
    # finds last '/' in string.
    subkey = event.symbol[:-event.symbol[::-1].find('/')-1]
    while subkey[-1] == '/': subkey = subkey[:-1]
    try: subdict = ip.user_ns["current_jobdict"][subkey] 
    except KeyError: raise IPython.ipapi.TryNext
    if hasattr(subdict, "children"): 
      if hasattr(subdict.children, "keys"):
        return [subkey + "/" + a + "/" for a in subdict.children.keys()]
    raise IPython.ipapi.TryNext
  else: return [a + "/" for a in ip.user_ns["current_jobdict"].children.keys()]

def current_jobname(self, arg):
  """ Returns current jobname. """
  ip = self.api
  if "current_jobdict" not in ip.user_ns: return
  print ip.user_ns["current_jobdict"].name
  return

def _main():
  import IPython.ipapi
  ip = IPython.ipapi.get()
  ip.expose_magic("explore", explore)
  ip.expose_magic("goto", goto)
  ip.expose_magic("listjobs", listjobs)
  ip.expose_magic("jobname", current_jobname)
  ip.set_hook('complete_command', goto_completer, re_key = '\s*%?goto')


_main()
