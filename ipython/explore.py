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
  from lada import interactive

  # options supported by all.
  parser = argparse.ArgumentParser(prog='%explore',
                     description='Opens a job-dictionary from file.')
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
  except SystemExit: return None
  else:
    if len(args.jobdict) == 0 and (args.type not in ["results", "errors", "all", "running"]):
      args.jobdict = args.type
      args.type = ""

  if     len(args.jobdict) == 0 \
     and (not args.is_file) \
     and (not args.is_expression) \
     and len(args.type) == 0 \
     and len(args.jobdict) == 0: 
    if interactive.jobdict is None:
      print "No current jobs."
    elif interactive.jobdict_path is None:
      print "Current position in job dictionary:", interactive.jobdict.name
    else:
      print "Current position in job dictionary:", interactive.jobdict.name
      print "Path to job dictionary: ", interactive.jobdict_path
    return

  options = ['', "errors", "results", "all"]
  if hasattr(self, "magic_qstat"): options.append("running")
  if args.type not in options: 
    print "Unknown TYPE argument {0}.\nTYPE can be one of {1}.".format(args.type, options)
    return

  # tries to open dictionary
  try: _explore_impl(self, args)
  except: return

  # now does special stuff if requested.
# if args.type == "errors": 
#   running_jobs = set(ip.magic("qstat").fields(-1)) if hasattr(self, "magic_qstat") else set([])
#   for name, job in interactive.jobdict.iteritems():
#     if interactive.jobdict_path is None: job.untag()
#     elif job.functional.Extract(join(dirname(interactive.jobdict_path),name)).success: job.tag()
#     elif name.replace("/", ".") in running_jobs: job.tag()
#     else: job.untag()

  if args.type == "results": 
    for name, job in interactive.jobdict.iteritems():
      if interactive.jobdict_path is None: job.tag()
      elif not job.functional.Extract(join(dirname(interactive.jobdict_path),name)).success: job.tag()
      else: job.untag()

# elif args.type == "running": 
#   running_jobs = set(ip.magic("qstat").fields(-1)) if hasattr(self, "magic_qstat") else set([])
#   for name, job in interactive.jobdict.iteritems():
#     if name.replace("/", ".") not in running_jobs: job.tag()
#     else: job.untag()

  elif args.type == "all": 
    for job in interactive.jobdict.itervalues(): job.untag()

def _explore_impl(self, args):
  """ Tries to open job-dictionary. """
  from os.path import abspath
  from copy import deepcopy
  from ..jobs import load, JobDict
  from ..jobs import JobParams, MassExtract as Collect
  from lada import interactive

  # case where we want to change the way the current dictionary is read.
  if len(args.jobdict) == 0:
    if interactive.jobdict is None:
      print "No job dictionary currently loaded.\n"\
            "Please use \"explore {0} path/to/jobict\".".format(args.type)
      interactive.__dict__.pop("current_jobdict_path", None)
      interactive.__dict__.pop("interactive.jobdict", None)
      return

    self.user_ns["interactive.jobdict"] = deepcopy(interactive.jobdict)
    if "collect" in self.user_ns: self.user_ns["collect"].uncache()
    interactive.__dict__.pop("_lada_subjob_iterator", None)
    interactive.__dict__.pop("_lada_subjob_iterated", None)
    return 

  # delete stuff from namespace.
  self.user_ns.pop("collect", None)
  self.user_ns.pop("jobparams", None)
  interactive.__dict__.pop("_lada_subjob_iterator", None)
  interactive.__dict__.pop("_lada_subjob_iterated", None)

  jobdict, new_path = None, None
  if args.is_file or not args.is_expression:
    try: jobdict = load(args.jobdict, timeout=60)
    except: jobdict = None
    else: new_path = abspath(args.jobdict)
  if jobdict is None and (args.is_expression or not args.is_file):
    try: jobdict = deepcopy(self.ev(args.jobdict))
    except: jobdict = None

  if not isinstance(jobdict, JobDict): jobdict = None

  if jobdict is None: # error
    print "Could not convert \"{0}\" to a job-dictionary.".format(args.jobdict) 
    return
    
  if args.concatenate: 
    interactive.jobdict.update(jobdict)
    jobdict = interactive.jobdict
    if new_path is None: new_path = interactive.jobdict_path
  elif args.add is not None:
    interactive.jobdict[args.add[1:-1]] = jobdict
    jobdict = interactive.jobdict
    if new_path is None: new_path = interactive.jobdict_path

  interactive.jobdict = jobdict
  self.user_ns["jobparams"] = JobParams()
  interactive.jobdict_path = new_path
  if new_path is not None: self.user_ns["collect"] = Collect(dynamic=True)

def completer(self, event): 
  """ Completer for explore. """ 
  from IPython import TryNext
  from IPython.core.completer import expand_user, compress_user
  from os.path import isdir
  from glob import iglob
  from ..jobs import JobDict
  from .. import jobdict_glob
  
  data = event.line.split()[1:]
  results, has_file, has_expr = [], False, False
  if "--file" in data: data.remove("--file"); has_file = True
  elif "--expression" in data: data.remove("--expression"); has_expr = True
  else: results = ["--file", "--expression"]
  if "--add" not in data: results.append('--add')
  if "--concatenate" not in data: results.append('--concatenate')
  
  results.extend(["errors", "results", "all"])
  if len(data) == 0: data = [''] 
  elif event.line[-1] == ' ': data.append('')
  if not has_file:
    results.extend( name for name, u in self.user_ns.iteritems()\
                    if isinstance(u, JobDict) and name[0] != '_' and name not in data)
  if not has_expr:
    relpath, tilde_expand, tilde_val = expand_user(data[-1])
    dirs = [f.replace('\\','/') + "/" for f in iglob(relpath+'*') if isdir(f)]
    dicts = [ f.replace('\\','/') for u in jobdict_glob for f in iglob(relpath+u)]
    if '.' in data[-1]:
      relpath, a, b = expand_user(data[-1][:data[-1].find('.')])
      dicts.extend([ f.replace('\\','/') for u in jobdict_glob for f in iglob(relpath+u)])
    dummy = [compress_user(p, tilde_expand, tilde_val) for p in dirs+dicts]
    results.extend(d for d in dummy if d not in data)
  return list(results)
