""" IPython explore function and completer. """
def explore(self, cmdl): 
  """ Starts exploration of a pickled job folder. 
  
      Usage: 
      The most standard form is to simply load a job folder. All other
      job-dictionary magic functions will then use it.
      
      >>> explore path/to/job_folder_pickle

      If you have created a job-folder directly (rather than save it to
      disk), you can also load it as

      >>> explore jobfolder_variable 

      In case of conflict between a pathname and a variable name, you can use
      the more explicit version.

      >>> explore --file jobfolder
      >>> explore --expression jobfolder

      You can load a dictionary and filter out successfull or unsuccessfull runs. 
      To explore errors only, use:
     
      >>> explore errors path/to/job_pickle

      To explore only successful results, use:

      >>> explore results path/to/job_pickle
  """

  import argparse
  from os.path import join, dirname
  from lada import interactive

  # options supported by all.
  parser = argparse.ArgumentParser(prog='%explore',
                     description='Opens a job-dictionary from file.')
  group = parser.add_mutually_exclusive_group()
  group.add_argument( '--file', action="store_true", dest="is_file",
                      help='JOBFOLDER is a path to a job-dictionary stored on disk.' )
  group.add_argument( '--expression', action="store_true", dest="is_expression",
                      help='JOBFOLDER is a python expression.' )
  parser.add_argument( 'type', metavar='TYPE', type=str, default="", nargs='?',
                       help="Optional. Specifies what kind of jobs will be explored. "\
                            "Can be one of results, errors, all, running. "\
                            "\"results\" are those jobs which have completed. "\
                            "\"errors\" are those jobs which are not \"running\" "\
                            "at the time of invokation and failed somehow. \"all\" means all jobs. "\
                            "By default, the dictionary is read as it was "\
                            "saved. The modified job-folder is not saved to "\
                            "disk.")
  parser.add_argument( 'jobfolder', metavar='JOBFOLDER', type=str, default="", nargs='?',
                       help='Job-dictionary variable or path to job folder saved to disk.')


  # parse arguments
  try: args = parser.parse_args(cmdl.split())
  except SystemExit: return None
  else:
    if len(args.jobfolder) == 0 and (args.type not in ["results", "errors", "all", "running"]):
      args.jobfolder = args.type
      args.type = ""

  if     len(args.jobfolder) == 0 \
     and (not args.is_file) \
     and (not args.is_expression) \
     and len(args.type) == 0 \
     and len(args.jobfolder) == 0: 
    if interactive.jobfolder is None:
      print "No current jobs."
    elif interactive.jobfolder_path is None:
      print "Current position in job folder:", interactive.jobfolder.name
    else:
      print "Current position in job folder:", interactive.jobfolder.name
      print "Path to job folder: ", interactive.jobfolder_path
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
  if args.type == "errors": 
    if interactive.jobfolder_path is None: 
      print "No known path/file for current job-folder.\n"\
            "Please save to file first."
      return
    running_jobs = set(self.magic("qstat").fields(-1)) if hasattr(self, "magic_qstat") else set([])
    for name, job in interactive.jobfolder.iteritems():
      if job.functional.Extract(join(dirname(interactive.jobfolder_path),name)).success: job.tag()
      elif name.replace("/", ".") in running_jobs: job.tag()
      else: job.untag()

  if args.type == "results": 
    if interactive.jobfolder_path is None: 
      print "No known path/file for current job-folder.\n"\
            "Please save to file first."
      return
    directory = dirname(interactive.jobfolder_path)
    for name, job in interactive.jobfolder.iteritems():
      if not job.functional.Extract(join(directory,name)).success: job.tag()
      else: job.untag()

# elif args.type == "running": 
#   running_jobs = set(ip.magic("qstat").fields(-1)) if hasattr(self, "magic_qstat") else set([])
#   for name, job in interactive.jobfolder.iteritems():
#     if name.replace("/", ".") not in running_jobs: job.tag()
#     else: job.untag()

  elif args.type == "all": 
    for job in interactive.jobfolder.itervalues(): job.untag()

def _explore_impl(self, args):
  """ Tries to open job-dictionary. """
  from os.path import abspath, isfile
  from ..jobs import load, JobFolder
  from ..jobs import JobParams, MassExtract as Collect
  from lada import interactive
  from lada.misc import LockFile, RelativePath

  # case where we want to change the way the current dictionary is read.
  if len(args.jobfolder) == 0:
    if interactive.jobfolder is None:
      print "No job folder currently loaded.\n"\
            "Please use \"explore {0} path/to/jobict\".".format(args.type)
      interactive.__dict__.pop("jobfolder", None)
      return

    if "collect" in self.user_ns: self.user_ns["collect"].uncache()
    interactive.__dict__.pop("_lada_subjob_iterator", None)
    interactive.__dict__.pop("_lada_subjob_iterated", None)
    return 

  # delete stuff from namespace.
  self.user_ns.pop("collect", None)
  self.user_ns.pop("jobparams", None)
  interactive.__dict__.pop("_lada_subjob_iterator", None)
  interactive.__dict__.pop("_lada_subjob_iterated", None)

  if args.is_file == False and args.is_expression == False \
     and isfile(RelativePath(args.jobfolder).path) \
     and isinstance(self.user_ns.get(args.jobfolder, None), JobFolder):
    print "The file {0} and the variable {1} both exist.\n"\
          "Please specify --file or --expression.\n"\
          .format(RelativePath(args.jobfolder).path, args.jobfolder)
    return
  jobfolder, new_path = None, None
  if args.is_file or not args.is_expression and isfile(RelativePath(args.jobfolder).path):
    try: jobfolder = load(args.jobfolder, timeout=6)
    except ImportError as e:
      print "ImportError: ", e
      return
    except Exception as e:
      print e
      if LockFile(args.jobfolder).is_locked:
        print "You may want to check for the existence of {0}."\
              .format(LockFile(args.jobfolder).lock_directory)
        print "If you are sure there are no jobs out there accessing {0},\n"\
              "you may want to delete that directory.".format(args.jobfolder)
      return
    else: new_path = abspath(args.jobfolder)
  if jobfolder is None and (args.is_expression or not args.is_file):
    jobfolder = self.user_ns.get(args.jobfolder, None)
    if not isinstance(jobfolder, JobFolder): 
      print "{0} is not a job-folder object.".format(args.jobfolder)
      return

  if jobfolder is None: # error
    print "Could not convert \"{0}\" to a job-dictionary.".format(args.jobfolder) 
    return
    
  interactive.jobfolder = jobfolder
  self.user_ns["jobparams"] = JobParams()
  interactive.jobfolder_path = new_path
  if new_path is not None: self.user_ns["collect"] = Collect(dynamic=True)

def completer(self, event): 
  """ Completer for explore. """ 
  from ..jobs import JobFolder
  from . import jobfolder_file_completer
  
  data = event.line.split()[1:]
  results, has_file, has_expr = [], False, False
  if "--file" in data: data.remove("--file"); has_file = True
  elif "--expression" in data: data.remove("--expression"); has_expr = True
  else: results = ["--file", "--expression"]
  
  results.extend(["errors", "results", "all"])
  if len(data) == 0: data = [''] 
  elif event.line[-1] == ' ': data.append('')
  if not has_file:
    results.extend( name for name, u in self.user_ns.iteritems()\
                    if isinstance(u, JobFolder) and name[0] != '_' and name not in data)
  if not has_expr: results.extend( jobfolder_file_completer(self, data))
  return list(results)
