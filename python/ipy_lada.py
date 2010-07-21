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
  
  args = [u for u in arg.split() if u not in ["in", "with"]]
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
      return result, abspath(splitpath(filename)[0]), splitpath(filename)[1]

    if is_a_variable: # checks the variable exists.
      var, dvar = None, ip.user_ns
      for name in  filename.split('.'):
        try: var, dvar = dvar[name], dvar[name].__dict__ 
        except:
          raise RuntimeError("Could not find variable name %s." % (filename))
      if not isinstance(var, jobs.JobDict): # checks the file exists.
        raise RuntimeError("Variable %s is not a JobDict object." % (filename))
      return var, None, None

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
      if "current_jobdict" not in ip.user_ns:
        raise RuntimeError("No job dictionary currently defined.\n"\
                           "Please load dictionary with \"explore\"." )
      if ip.user_ns["current_jobdict_directory"] == None: 
        raise RuntimeError("Directory of current job dictionary is not defined.\n"\
                           "Cannot check for success.\n"
                           "Please reload dictionary with explore, \n"
                           "or set current_jobdict_directory by hand.")
      ip.user_ns["current_jobdict_filename"] = None
      ip.user_ns["current_jobdict"] = deepcopy(ip.user_ns["current_jobdict"])
      rootdir = ip.user_ns["current_jobdict_directory"]
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
    elif len(args) == 1:
      result = get_dict(args[0])
      if result != None:
        ip.user_ns["current_jobdict"] = result[0]
        ip.user_ns["current_jobdict_directory"] = result[1]
        ip.user_ns["current_jobdict_filename"] = result[2]
    elif len(args) == 2 and args[0] == "file":
      result = get_dict(args[0], True)
      if result != None:
        ip.user_ns["current_jobdict"] = result[0]
        ip.user_ns["current_jobdict_directory"] = result[1]
        ip.user_ns["current_jobdict_filename"] = result[2]
    elif len(args) == 2 and args[0] == "JobDict":
      result = get_dict(args[0], False)
      if result != None:
        ip.user_ns["current_jobdict"] = result[0]
        ip.user_ns["current_jobdict_directory"] = result[1]
        ip.user_ns["current_jobdict_filename"] = result[2]
    else: raise RuntimeError("Calling explore with arguments %s is invalid." % (arg))

  try: _impl()
  except RuntimeError as e: print e
  else:
    # remove any prior iterator stuff.
    ip.user_ns.pop("_lada_subjob_iterator", None)
    ip.user_ns.pop("_lada_subjob_iterated", None)

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
    if arg == "all": 
      for job, d in current.root.walk_through():
        if job.is_tagged: continue
        print job.name
      return
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
  if "current_jobdict" not in ip.user_ns:
    print "No current job-dictionary. Please load one with \"explore\"."
    return

  args = event.split()
  if len(args) > 1: 
    print "Invalid argument %s." % (event)
  elif len(args) == 0:
    if "_lada_subjob_iterator" in ip.user_ns: 
      iterator = ip.user_ns["_lada_subjob_iterator"]
    else:
      iterator = ip.user_ns["current_jobdict"].root.walk_through()
    while True:
      try: job, d = iterator.next()
      except StopIteration: 
        print "Reached end of job list."
        return 
      if job.is_tagged: continue
      break
    ip.user_ns["_lada_subjob_iterator"] = iterator
    if "_lada_subjob_iterated" not in ip.user_ns: ip.user_ns["_lada_subjob_iterated"] = []
    ip.user_ns["_lada_subjob_iterated"].append(ip.user_ns["current_jobdict"].name)
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
      print "In job ", ip.user_ns["current_jobdict"].name


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
  else:
    result = [a + "/" for a in ip.user_ns["current_jobdict"].children.keys()]
    result.extend(["/", "next", "reset"])
    if ip.user_ns["current_jobdict"].parent != None: result.append("../")
    if "_lada_subjob_iterated" in ip.user_ns:
      if len(ip.user_ns["_lada_subjob_iterated"]): result.append("previous")
    return result

def showme(self, event):
  """ Edits functional and/or structure. """
  from os import remove
  from tempfile import NamedTemporaryFile
  from lada.opt import read_input
  from lada.crystal import write_poscar, read_poscar
  ip = self.api
  if "current_jobdict" not in ip.user_ns: 
    print "No job dictionary. Please use \"explore\" magic function."
    return
  args = [u for u in event.split() if u not in ["in"]]
  if len(args) == 0:
    print "What should I show you?"
    return
  if len(args) == 2:
    old = ip.user_ns["current_jobdict"].name 
    try:
      goto(self, args[0])
      showme(self, args[1])
    finally: goto(self, old)
  elif len(args) == 1:
    arg = args[0]
    job = ip.user_ns["current_jobdict"]
    try: 
      suffix = ".py" if arg.lower() == "functional" else None
      with NamedTemporaryFile("w", delete = False, suffix=suffix) as file:
        filename = file.name
        # editing INCAR.
        if arg.lower() == "functional":
          if job.functional == None: # case where no job is defined.
            file.write("# There are currently no actual jobs defined here.\n"\
                       "functional = None\njobparams={}\n")
          else: # case where a functional exists.
            file.write(repr(job.functional))
            file.write("\n\n# Parameters in the functional above **will** be \n"\
                       "# overwritten by the following corresponding parameters.\n")
            file.write("jobparams = {}\n")
            for key, value in job.jobparams.items():
              if key == "args": continue
              if key == "functional": continue
              file.write("jobparams[\"%s\"] = %s\n" % (key, repr(value)))
        # editing POSCAR.
        elif arg.lower() == "structure":
          # always as vasp5. Makes it easier.
          if job.args[0] != None: write_poscar(job.args[0], file, True)
        # Error!
        else: 
          print "%s is not a valid argument to showme." % (event)
          return

      # lets user edit stuff.
      ip.magic("%%edit -x %s" % (filename))

      # change jobparameters.
      if arg.lower() == "functional":
        input = read_input(filename)
        input.jobparams["args"] = job.args
        input.jobparams["functional"] = input.functional
        job.jobparams = input.jobparams
      elif arg.lower() == "structure": 
        job.args[0] = read_poscar(path=filename)
    finally: remove(filename)

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
