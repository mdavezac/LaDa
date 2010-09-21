""" IPython explore function and completer. """
def explore(self, cmdl): 
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

      >>> explore --file jobdict
      >>> explore --variable jobdict

      You can load a dictionary and filter out successfull or unsuccessfull runs. 
      To explore errors only, use:
     
      >>> explore errors in path/to/job_pickle

      To explore only successful results, use:

      >>> explore results in path/to/job_pickle
  """

  import argparse
  from os.path import join, dirname
  from . import _get_current_job_params

  ip = self.api
  # removes current error if any.
  ip.user_ns.pop("_lada_error", None)
  
  # options supported by all.
  parser = argparse.ArgumentParser(prog='%explore')
  group = parser.add_mutually_exclusive_group()
  group.add_argument( '--file', action="store_true", dest="is_file",
                      help='JOBDICT is a path to a job-dictionary stored on disk.' )
  group.add_argument( '--expression', action="store_true", dest="is_expression",
                      help='JOBDICT is a python expression.' )
  group.add_argument( '--concatenate', action="store_true", dest="concatenate",
                      help='Update/overwrites current jobdictionary with newly explored one.' )
  group.add_argument( '--add', metavar="NAME", type=str, dest="add",
                      help='Adds newly explored jobdictionary to current'\
                           'dictionary in specified location..' )
  parser.add_argument( 'type', metavar='TYPE', type=str, default="", nargs='?',
                       help="Optional. Specifies what kind of jobs will be explored. "\
                            "Can be one of results, errors, all, running. "\
                            "\"results\" are those jobs which have completed. "\
                            "\"errors\" are those jobs which are not \"running\" "\
                            "at the time of invokation and failed somehow. \"all\" means all jobs. "\
                            "By default, the dictionary is read as it was saved. ")
  parser.add_argument( 'jobdict', metavar='JOBDICT', type=str, default="", nargs='?',
                       help='Job-dictionary variable or path to job dictionary saved to disk.')


  # parse arguments
  try: args = parser.parse_args(cmdl.split())
  except SystemExit as e: return None
  else:
    if len(args.jobdict) == 0 and (args.type not in ["results", "errors", "all", "running"]):
      args.jobdict = args.type
      args.type = ""

  if     len(args.jobdict) == 0 \
     and (not args.is_file) \
     and (not args.is_expression) \
     and len(args.type) == 0 \
     and len(args.jobdict) == 0: 
    current, path = _get_current_job_params(self, 0)
    if current == None:
      print "No current jobs."
    elif path == None:
      print "Current position in job dictionary:", current.name
    else:
      print "Current position in job dictionary:", current.name
      print "Path to job dictionary: ", path
    return

  options = ['', "errors", "results", "all"]
  if hasattr(self, "magic_qstat"): options.append("running")
  if args.type not in options: 
    ip.user_ns["_lada_error"] = "Unknown TYPE argument {0}.\n"\
                                "TYPE can be one of {1}.".format(args.type, options)
    print ip.user_ns["_lada_error"]
    return

  # tries to open dictionary
  try: _explore(self, args)
  except: return

  # now does special stuff if requested.
  current, path = _get_current_job_params(self, 0)
  if args.type == "errors": 
    running_jobs = set(ip.magic("qstat").fields(-1)) if hasattr(self, "magic_qstat") else set([])
    for job, name in current.walk_through():
      if path == None: job.untag()
      elif job.functional.Extract(join(dirname(path),name)).success: job.tag()
      elif name.replace("/", ".") in running_jobs: job.tag()
      else: job.untag()

  elif args.type == "results": 
    for job, name in current.walk_through():
      if path == None: job.tag()
      elif not job.functional.Extract(join(dirname(path),name)).success: job.tag()
      else: job.untag()

  elif args.type == "running": 
    running_jobs = set(ip.magic("qstat").fields(-1)) if hasattr(self, "magic_qstat") else set([])
    for job, name in current.walk_through():
      if name.replace("/", ".") not in running_jobs: job.tag()
      else: job.untag()

  elif args.type == "all": 
    for job, name in current.walk_through(): job.untag()

def explore_completer(self, event): 
  """ Completer for explore. """ 
  import IPython

  ip = self.api

  options = set(["errors", "results", "all", "--file", "--expression"])
  if hasattr(self, "magic_qstat"): options.add("running")
  data = event.line.split()[1:]
  if len(set(data) - options) > 1: raise IPython.ipapi.TryNext

  results, has_file, has_expr = [], False, False
  if "--file" in data: data.remove("--file"); has_file = True
  elif "--expression" in data: data.remove("--expression"); has_expr = True
  else: results = ["--file", "--expression"]
  if "--add" not in data: results.append('--add')
  if "--concatenate" not in data: results.append('--concatenate')

  results.extend(["errors", "results", "all"])
  if hasattr(self, "magic_qstat"): results.append("running")
  if len(data) == 0: 
    results.extend(_glob_job_pickles(ip, ""))
    return results

  if len(event.symbol) == 0: data.append("")
  
  if not has_file:
    results.extend(_glob_ipy_user_ns(ip, data[-1]) )
  if not has_expr:
    results.extend(_glob_job_pickles(ip, data[-1]) )
  return results


def _glob_job_pickles(ip, arg):
  """ Returns list of candidate files and directory. """
  from os.path import join

  # nothing there yet.
  try: 
    result = [ u + "/" for u in ip.magic("mglob dir:{0}*".format(arg)) ]
    s = "mglob \"cont:JobDict pickle\" {0}*".format(arg)
    result.extend([u for u in ip.magic(s)])
    return result
  except: print "Error in _explore._glob_job_pickles."

  return []


def _glob_ipy_user_ns(ip, arg):
  """ Returns list of candidate python variables. """
  from ..jobs import JobDict
  
  if len(arg) == 0:
    return [u for u in ip.user_ns if isinstance(ip.user_ns[u], JobDict)]

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


def _explore(self, args):
  """ Tries to open job-dictionary. """
  from os.path import exists, abspath, dirname
  from copy import deepcopy
  from ..jobs import MassExtract, load, JobDict
  from ._collect import Collect
  from . import _get_current_job_params

  ip = self.api

  # case where we want to change the way the current dictionary is read.
  if len(args.jobdict) == 0:
    c = ip.user_ns.pop("current_jobdict", None)
    if c == None:
      ip.user_ns["_lada_error"] = "No job dictionary currently loaded.\n"\
                                  "Please use \"explore {0} path/to/jobict\".".format(args.type)
      ip.user_ns.pop("current_jobdict_path", None)
      print ip.user_ns["_lada_error"] 
      raise RuntimeError(ip.user_ns["_lada_error"])

    ip.user_ns["current_jobdict"] = deepcopy(c)
    if "collect" in ip.user_ns: ip.user_ns["collect"].uncache()
    ip.user_ns.pop("_lada_subjob_iterator", None)
    ip.user_ns.pop("_lada_subjob_iterated", None)
    return 

  # delete stuff from namespace.
  current, path = _get_current_job_params(self, 0)
  ip.user_ns.pop("collect", None)
  ip.user_ns.pop("_lada_subjob_iterator", None)
  ip.user_ns.pop("_lada_subjob_iterated", None)

  jobdict, new_path = None, None
  if args.is_file or not args.is_expression:
    try: jobdict = load(args.jobdict)
    except: jobdict = None
    else: new_path = abspath(args.jobdict)
  if jobdict == None and (args.is_expression or not args.is_file):
    try: jobdict = deepcopy(ip.ev(args.jobdict))
    except: jobdict = None

  if not isinstance(jobdict, JobDict): jobdict = None

  if jobdict == None: # error
    ip.user_ns["_lada_error"] = \
       "Could not convert \"{0}\" to a job-dictionary.".format(args.jobdict) 
    print ip.user_ns["_lada_error"] 
    raise RuntimeError(ip.user_ns["_lada_error"])
    
  if args.concatenate: 
    current.update(jobdict)
    jobdict = current
    if new_path == None: new_path = path
  elif args.add != None:
    current[args.add[1:-1]] = jobdict
    jobdict = current
    if new_path == None: new_path = path

  ip.user_ns["current_jobdict"] = jobdict
  if new_path != None: ip.user_ns["current_jobdict_path"] = new_path
  else: ip.user_ns.pop("current_jobdict_path", None)
  if new_path != None: ip.user_ns["collect"] = Collect()
