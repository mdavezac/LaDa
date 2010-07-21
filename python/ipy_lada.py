""" IPython functions and data. """
from . import jobs

def _get_current_job_params(self, verbose=0):
  """ Returns a tuple with current job, filename, directory. """

  ip = self.api
  if "current_jobdict" not in ip.user_ns: 
    if verbose > 0:
      print "No job dictionary. Please use \"explore\" magic function."
    return None, None
  current = ip.user_ns["current_jobdict"]
  if "current_jobdict_path" not in ip.user_ns: 
    if verbose > 1:
      print "No filepath for current job dictionary.\n"\
            "Please set current_jobdict_path."
    return None, None, None
  path = ip.user_ns["current_jobdict_path"]
  return current, path



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
  from os.path import exists, split as splitpath, abspath
  from cPickle import load
  from .opt import open_exclusive
  ip = self.api
  current, path = _get_current_job_params(self, 0)
  
  args = [u for u in arg.split() if u not in ["in", "with"]]
  if len(args) == 0:
    if current == None:
      print "No current jobs."
    elif path == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Path to job dictionary: ", path

  if isinstance(arg, jobs.JobDict):
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
      with open_exclusive(filename, "r") as file: result = load(file)
      return result, abspath(filename)

    if is_a_variable: # checks the variable exists.
      var, dvar = None, ip.user_ns
      for name in  filename.split('.'):
        try: var, dvar = dvar[name], dvar[name].__dict__ 
        except:
          raise RuntimeError("Could not find variable name %s." % (filename))
      if not isinstance(var, jobs.JobDict): # checks the file exists.
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
  except RuntimeError as e: print e
  else:
    # remove any prior iterator stuff.
    ip.user_ns.pop("_lada_subjob_iterator", None)
    ip.user_ns.pop("_lada_subjob_iterated", None)

def goto(self, arg):
  """ Moves current dictionary position and working directory (if appropriate). """
  from os import chdir
  from os.path import exists, join, split as splitpath
  ip = self.api
  current, path = _get_current_job_params(self, 1)
    
  if current == None: 
    print "No current jobs."
    return
  if len(arg.split()) == 0:
    if path == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Filename of jobdictionary: ", path
    return
  args = arg.split()
  if len(args) > 1:
    print "Invalid argument to goto %s." % (arg)
    return

  # if no argument, then print current job data.
  if len(args) == 0: 
    explore(self, "")
    return

  # cases to send to iterate.
  if args[0] == "next":       return iterate(self, "")
  elif args[0] == "previous": return iterate(self, "previous")
  elif args[0] == "reset":    return iterate(self, "reset")

  # case for which precise location is given.
  try: result = current[args[0]] 
  except KeyError as e: 
    print e
    return 

  if not (hasattr(result, "parent") and hasattr(result, "children")\
     and hasattr(result, "jobparams")):
    print "%s is a job parameter, not a directory."
    return
  ip.user_ns["current_jobdict"] = result
  current = ip.user_ns["current_jobdict"]
  if path == None: return
  dir = join(splitpath(path)[0], current.name[1:]) 
  if exists(dir): chdir(dir)
  return

def listjobs(self, arg):
  """ Lists subjobs. """
  ip = self.api
  current, path = _get_current_job_params(self, 1)
  if current == None: return
  if len(arg) != 0:
    if arg == "all": 
      for job, d in current.root.walk_through():
        if job.is_tagged: continue
        print job.name
      return
    try: subdict = current[arg] 
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
    if current.children[j].is_tagged: continue
    if len(string+j) > 60:
      if len(lines) != 0: lines += "\n" + string
      else: lines = string
      string = ""
    string += j + " "

  if len(lines) == 0: print string
  else: print lines + "\n" + string

def iterate(self, event):
  """ Goes to next (untagged) job. """
  ip = self.api
  current, path = _get_current_job_params(self, 1)
  if current == None: return

  args = event.split()
  if len(args) > 1: 
    print "Invalid argument %s." % (event)
  elif len(args) == 0:
    if "_lada_subjob_iterator" in ip.user_ns: 
      iterator = ip.user_ns["_lada_subjob_iterator"]
    else:
      iterator = current.root.walk_through()
    while True:
      try: job, d = iterator.next()
      except StopIteration: 
        print "Reached end of job list."
        return 
      if job.is_tagged: continue
      break
    ip.user_ns["_lada_subjob_iterator"] = iterator
    if "_lada_subjob_iterated" not in ip.user_ns: ip.user_ns["_lada_subjob_iterated"] = []
    ip.user_ns["_lada_subjob_iterated"].append(current.name)
    goto(self, job.name)
    print "In job ", ip.user_ns["current_jobdict"].name
  elif args[0] == "reset" or args[0] == "restart":
    # remove any prior iterator stuff.
    ip.user_ns.pop("_lada_subjob_iterator", None)
    ip.user_ns.pop("_lada_subjob_iterated", None)
    print "In job ", ip.user_ns["current_jobdict"].name
  elif args[0] == "back" or args[0] == "previous":
    if "_lada_subjob_iterated" not in ip.user_ns: 
      print "No previous job to go to. "
    else:
      goto(self, ip.user_ns["_lada_subjob_iterated"].pop(-1))
      if len(ip.user_ns["_lada_subjob_iterated"]) == 0:
        del ip.user_ns["_lada_subjob_iterated"]
      print "In job ", ip.user_ns["current_jobdict"].name


def goto_completer(self, event):
  import IPython
  ip = self.api
  current, path = _get_current_job_params(self, 0)
  if current == None: raise IPython.ipapi.TryNext

  elif '/' in event.symbol:
    # finds last '/' in string.
    subkey = event.symbol[:-event.symbol[::-1].find('/')-1]
    while subkey[-1] == '/': subkey = subkey[:-1]
    try: subdict = current[subkey] 
    except KeyError: raise IPython.ipapi.TryNext
    if hasattr(subdict, "children"): 
      if hasattr(subdict.children, "keys"):
        return [subkey + "/" + a + "/" for a in subdict.children.keys()]
    raise IPython.ipapi.TryNext
  else:
    result = [a + "/" for a in current.children.keys()]
    result.extend(["/", "next", "reset"])
    if current.parent != None: result.append("../")
    if "_lada_subjob_iterated" in ip.user_ns:
      if len(ip.user_ns["_lada_subjob_iterated"]): result.append("previous")
    return result

def showme(self, event):
  """ Edits functional and/or structure. """
  from os import remove, stat
  from tempfile import NamedTemporaryFile
  from lada.opt import read_input
  from lada.crystal import write_poscar, read_poscar, Structure
  ip = self.api
  # gets dictionary, path.
  current, path = _get_current_job_params(self, 1)
  if current == None: return
  # splits argumetns, removes decorative keywords.
  args = [u for u in event.split() if u not in ["in"]]
  # nothing to do.
  if len(args) == 0: 
    print "What should I show you?"
    return
  # showme *whatever* in *where*
  if len(args) == 2:
    old = current.name 
    try:
      goto(self, args[1]) # goto *where*
      showme(self, args[0]) # show *whatever*
    finally: goto(self, old) # go back.
  # showme *whatever*
  elif len(args) == 1:
    arg = args[0]
    filename = None
    try: # try/finally section will removed namedtemporaryfile.
      # want .py suffix to get syntax highlighting in editors.
      suffix = ".py" if arg.lower() == "functional" else None
      with NamedTemporaryFile("w", delete = False, suffix=suffix) as file:
        filename = file.name
        # editing INCAR.
        if arg.lower() == "functional":
          if current.functional == None: # case where no job is defined.
            file.write("# There are currently no actual jobs defined here.\n"\
                       "functional = None\njobparams={}\n")
          else: # case where a functional exists.
            file.write(repr(current.functional))
            file.write("\n\n# Parameters in the functional above **will** be \n"\
                       "# overwritten by the following corresponding parameters.\n")
            file.write("jobparams = {}\n")
            for key, value in current.jobparams.items():
              if key == "args": continue
              if key == "functional": continue
              file.write("jobparams[\"%s\"] = %s\n" % (key, repr(value)))
        # editing POSCAR.
        elif arg.lower() == "structure":
          # always as vasp5. Makes it easier.
          structure = Structure() if current.args[0] == None else current.args[0] 
          write_poscar(current.args[0], file, True)
        # Error!
        else: 
          print "%s is not a valid argument to showme." % (event)
          return

      # lets user edit stuff.
      time0 = stat(filename)[-2]
      ip.magic("%%edit -x %s" % (filename))
      if stat(filename)[-2] == time0: return

      # change jobparameters.
      if arg.lower() == "functional":
        input = read_input(filename)
        input.jobparams["args"] = job.args
        input.jobparams["functional"] = input.functional
        current.jobparams = input.jobparams
      elif arg.lower() == "structure": 
        current.args[0] = read_poscar(path=filename)
    finally:
      if filename != None:
        try: remove(filename)
        except: pass


def saveto(self, event):
  """ Saves current job to current filename and directory. """
  from os.path import exists, abspath, isfile
  from lada import jobs
  ip = self.api
  # gets dictionary, path.
  current, path = _get_current_job_params(self, 1)
  if current == None: return
  args = [u for u in event.split() ]
  if len(args) == 0: 
    if path == None: 
      print "No current job-dictionary path.\n"\
            "Please specify on input, eg"\
            ">saveto this/path/filename"
      return
    if exists(path): 
      if not isfile(path): 
        print "%s is not a file." % (path)
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File already exists. Overwrite? [y/n] ")
      if a == 'n': return
    jobs.save(current, path) 
  elif len(args) == 1:
    if exists(args[0]): 
      if not isfile(path): 
        print "%s is not a file." % (path)
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File already exists. Overwrite? [y/n] ")
      if a == 'n': return
    jobs.save(current, args[0]) 
    ip.user_ns["current_jobdict_path"] = abspath(args[0])
  else: print "Invalid call to saveto."


def showme_completer(self, event):
  import IPython
  ip = self.api
  if "current_jobdict" not in ip.user_ns: raise IPython.ipapi.TryNext
  return ["structure", "functional"]

def current_jobname(self, arg):
  """ Returns current jobname. """
  ip = self.api
  if "current_jobdict" not in ip.user_ns: return
  print ip.user_ns["current_jobdict"].name
  return

def qstat(self, arg):
  """ squeue --user=`whoami` -o "%7i %.3C %3t  --   %50j" """
  from subprocess import Popen, PIPE
  from IPython.genutils import SList

  ip = self.api
  # finds user name.
  whoami = Popen(["whoami"], stdout=PIPE).stdout.readline()[:-1]
  squeue = Popen(["squeue", "--user=" + whoami, "-o", "\"%7i %.3C %3t    %j\""],
                 stdout=PIPE)
  result = squeue.stdout.read().rstrip().split('\n')
  result = SList([u[1:-1] for u in result])
  return result.grep(str(arg[1:-1]))

def cancel_jobs(self, arg):
  """ Cancel jobs which grep for whatever is in arg.
  
      For instance, the following cancels all jobs with "anti-ferro" in their
      name.
      >>> %cancel_jobs "anti-ferro"
  """
  from subprocess import Popen, PIPE
  
  arg = str(arg[1:-1])
  if len(arg) == 0: 
    print "cancel_job Requires an argument."
    print "Please use please_cancel_all_jobs to cancel all jobs."
    return
  result = qstat(self, arg)
  for u, name in zip(result.fields(0), result.fields(-1)):
    print "cancelling %s." % (name)
  a = ''
  while a not in ['n', 'y']:
    a = raw_input("Are you sure you want to cancel all jobs? [y/n] ")
  if a == 'n': return
  for u, name in zip(result.fields(0), result.fields(-1)):
    ip.system("scancel %i" % (int(u)))

def please_cancel_all_jobs(self, arg):
  """ Cancel all jobs. """
  from subprocess import Popen, PIPE
  
  a = ''
  while a not in ['n', 'y']: a = raw_input("Are you sure you want to cancel all jobs? [y/n] ")
  if a == 'n': return
  result = qstat(self, None)
  for u in result.field(0):
    ip.system("scancel %i" % (int(u)))


def _main():
  import lada
  from os import environ
  import IPython.ipapi

  ip = IPython.ipapi.get()
  ip.expose_magic("explore", explore)
  ip.expose_magic("goto", goto)
  ip.expose_magic("listjobs", listjobs)
  ip.expose_magic("jobname", current_jobname)
  ip.expose_magic("iterate", iterate)
  ip.expose_magic("showme", showme)
  ip.expose_magic("savejobs", saveto)
  ip.set_hook('complete_command', goto_completer, re_key = '\s*%?goto')
  ip.set_hook('complete_command', showme_completer, re_key = '\s*%?showme')
  if "SNLCLUSTER" in environ:
    if environ["SNLCLUSTER"] in ["redrock"]:
      ip.expose_magic("qstat", qstat)
      ip.expose_magic("cancel_jobs", cancel_jobs)
      ip.expose_magic("please_cancel_all_jobs", please_cancel_all_jobs)

  for key in lada.__dict__:
    if key[0] == '_': continue
    if key == "jobs": ip.ex("from lada import jobs as ladajobs")
    else: ip.ex("from lada import " + key)


_main()
