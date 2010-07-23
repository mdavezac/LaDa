""" IPython functions and data. """
from contextlib  import contextmanager

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


def goto(self, arg):
  """ Moves current dictionary position and working directory (if appropriate). """
  from os import chdir
  from os.path import exists, join, split as splitpath
  ip = self.api
  current, path = _get_current_job_params(self, 1)
  ip.user_ns.pop("_lada_error", None)
    
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
    ip.user_ns["_lada_error"] = "Invalid argument to goto %s." % (arg)
    print ip.user_ns["_lada_error"]
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
    ip.user_ns["_lada_error"] = e
    return 

  if not (hasattr(result, "parent") and hasattr(result, "children")\
     and hasattr(result, "jobparams")):
    ip.user_ns["_lada_error"] = "%s is a job parameter, not a directory."
    print ip.user_ns["_lada_error"] 
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
  ip.user_ns.pop("_lada_error", None)
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
  ip.user_ns.pop("_lada_error", None)
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
          structure = Structure() if current.structure == None else current.structure 
          write_poscar(current.structure, file, True)
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
        current.structure = read_poscar(path=filename)
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
  ip.user_ns.pop("_lada_error", None)
  if current == None:
    ip.user_ns["_lada_error"] = "No job-dictionary to save."
    print ip.user_ns["_lada_error"] 
    return
  args = [u for u in event.split() ]
  if len(args) == 0: 
    if path == None: 
      ip.user_ns["_lada_error"] = "No current job-dictionary path.\n"\
                                  "Please specify on input, eg"\
                                  ">saveto this/path/filename"
      print ip.user_ns["_lada_error"] 
      return
    if exists(path): 
      if not isfile(path): 
        ip.user_ns["_lada_error"] = "%s is not a file." % (path)
        print ip.user_ns["_lada_error"] 
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File already exists. Overwrite? [y/n] ")
      if a == 'n': return
    jobs.save(current, path, overwrite=True) 
  elif len(args) == 1:
    if exists(args[0]): 
      if not isfile(path): 
        ip.user_ns["_lada_error"] = "%s is not a file." % (path)
        print ip.user_ns["_lada_error"] 
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File already exists. Overwrite? [y/n] ")
      if a == 'n': return
    jobs.save(current, args[0], overwrite=True) 
    ip.user_ns["current_jobdict_path"] = abspath(args[0])
  else:
    ip.user_ns["_lada_error"] = "Invalid call to saveto."
    print ip.user_ns["_lada_error"] 


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

@contextmanager
def save_state(self):
  from copy import deepcopy
  ip = self.api
  current, path = _get_current_job_params(self, 0)
  if current != None: current.copy()
  iterated, iterator = None, None
  if _lada_subjob_iterator in ip.user_ns:
    iterated = deepcopy(ip.user_ns["_lada_subjob_iterator"])
  if _lada_subjob_iterator in ip.user_ns:
    iterator = deepcopy(ip.user_ns["_lada_subjob_iterator"])
  yield current, path, iterated, iterator

  if current == None: ip.user_ns.pop("current_jobdict", None)
  else: ip.user_ns["current_jobdict"] = current
  if path == None: ip.user_ns.pop("current_jobdict_path", None)
  else: ip.user_ns["current_jobdict_path"] = path
  if iterator == None: ip.user_ns.pop("_lada_subjob_iterator", None)
  else: ip.user_ns["_lada_subjob_iterator"] = iterator
  if iterated == None: ip.user_ns.pop("_lada_subjob_iterated", None)
  else: ip.user_ns["_lada_subjob_iterated"] = iterated



def launch_jobs(self, event):
  """ Launch pbs jobs according to given recipe. """
  import re
  from os.path import split as splitpath, join, exists
  from lada.opt.changedir import Changedir
  ip = self.api
  current, path = _get_current_job_params(self, 0)
  ip.user_ns.pop("_lada_error", None)

  # parse input phrase 
  scattered = re.search("scattered", event)
  mppalloc_re = re.search("\s+(?:with)?\s+mppalloc\s*=?\s*(\S+)", event);
  pooled = re.search("(?:in)?\s*(\d+)?\s+pools\s+(?:of)?\s+(\d+)\s+(?:processors)?", event)
  filepath = re.search("(?:as)?\s+filepath\s+(\S+)", event)
  jobdict = re.search("^\s*(file\s|JobDict\s)\s*(\S+)")
  andlaunch = re.search("(?:and\s)?\s*launch\s*$")
  if scattered != None and pooled != None: 
    print "Invalid argument.\n"\
          "Choose either pooled or scattered script parallelization method."
    return 
  if pooled != None and mppalloc != None:
    print "Invalid argument.\n"\
          "mppalloc argument is meaningless in the context of pooled pbs parallelization."
    return
  if pooled == None and scattered == None:
    print "invalid argument.\n"
    return
  if jobdict != None: # make sure we are not catching scattered or pooled directive.
    if jobdict.group(2) == "scattered": jobdict == None
    elif jobdict.group(2) == "pooled": jobdict == None
    else:
      try: int(jobdict.group(2))
      except: pass
      else: jobdict = None

  # gets filepath
  if filepath != None: filepath = filepath.match(1)

  # save current state.
  with save_state() as state:
    # opens new dictionary if requested.
    if jobdict != None:
      explore(self, jobdict.group(0))
      if "_lada_error" in ip.user_ns:
        print "Encountered error. Will not launch jobs."
        return
      current, path = _get_current_job_params(self, 0)
    # saves to new path if requested.
    if filepath != None: save(self, filepath)
    elif jobdict == None: save(self, "")
    if "_lada_error" in ip.user_ns:
      print "Encountered error. Will not launch jobs."
      return
    # required scattered parallelization (one pbs script per job).
    if scattered != None: 
      # creates mppalloc function.
      def mppalloc(job): 
        """ Returns number of processes for this job. """
        N = len(job.structure.atoms) # number of atoms.
        if N % 2 == 1: N -= 1
        return N  
      # gets used defined mppalloc if requested.
      if mppalloc_re != None:
        try: mppalloc = int(mppalloc_re.group(1)) 
        except ValueError:
          if mppalloc_re.group(1) not in ip.user_ns:
            print "Could not find processor allocation function %s." % (mppalloc.group(1))
            return
          mppalloc = ip.user_ns[mppalloc_re.group(1)]
      if path == None: 
        print "Could not determine which directory to use for output."
        return

      # where are we? important for default template.
      which = "SNLCLUSTER" in environ
      if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
      template = default_slurm if which else default_pbs
      # gets python script to launch in pbs.
      pyscript = __file__.replace(splitpath(__file__)[1], "runone.py")
      # creates directory.
      with Changedir(path + ".pbs") as pwd: pass 
      # creates pbs scripts.
      pbsscripts = []
      for i, (job, name) in enumerate(current.walk_through()):
        if job.is_tagged: continue
        mppwidth = mppalloc(job) if hasattr(mppalloc, "__call__") else mppalloc
        name = name.replace("/", ".")
        pbsscripts.append( abspath(join(directory, name + ".pbs")) )
        with open(results[-1], "w") as file: 
          jobs.template( file, outdir=splitpath(path[0]), jobid=i, mppwidth=mppwidth, name=name,\
                         pickle = path, pyscript=pyscript )
      # return if not launching
      if andlaunch == None: return result
      # otherwise, launch.
      for script in pbsscripts:
        if which: ip.ex("qsub %s" % script)
        else: ip.ex("sbatch %s" % script)

    # required pooled parallelization (few pbs scripts, many jobs)
    else: 
      print "*pooled* not yet implemented."
      return


def fakerun(self, event):
  """ Creates job directory tree and input files without computing. """
  from os.path import path, split as splitpath

  current, path = _get_current_job_params(self, 0)
  ip.user_ns.pop("_lada_error", None)
  if len(event.split()) > 1: 
    print "fakerun does not take an argument."
    return
  elif len(event.split()) == 1: directory = event.split()[0]
  else: directory = splitpath(path)[0]

  if exists(directory) and not isdir(directory):
    print "%s exists and is not a directory." % (directory)
    return 
  elif exists(directory):
    a = ''
    while a not in ['n', 'y']:
      a = raw_input("%s exists. Some input files could be overwritten. Continue? [y/n]")
    if a == 'n': return
  for job, dirname in jobdict.walk_through(path):
    if not job.is_tagged: job.compute(outdir=dirname, norun=True)


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
  from os import environ
  import IPython.ipapi
  import lada
  from ._goto import goto, iterate, goto_completer
  from ._explore import explore, explore_completer
  from ._showme import showme, showme_completer

  ip = IPython.ipapi.get()
  ip.expose_magic("explore", explore)
  ip.expose_magic("goto", goto)
  ip.expose_magic("listjobs", listjobs)
  ip.expose_magic("jobname", current_jobname)
  ip.expose_magic("iterate", iterate)
  ip.expose_magic("showme", showme)
  ip.expose_magic("savejobs", saveto)
  ip.expose_magic("fakerun", fakerun)
  ip.set_hook('complete_command', goto_completer, re_key = '\s*%?goto')
  ip.set_hook('complete_command', showme_completer, re_key = '\s*%?showme')
  ip.set_hook('complete_command', explore_completer, re_key = '\s*%?explore')
  if "SNLCLUSTER" in environ:
    if environ["SNLCLUSTER"] in ["redrock"]:
      ip.expose_magic("qstat", qstat)
      ip.expose_magic("cancel_jobs", cancel_jobs)
      ip.expose_magic("please_cancel_all_jobs", please_cancel_all_jobs)

  for key in lada.__dict__:
    if key[0] == '_': continue
    if key == "ipython": continue
    if key == "jobs": ip.ex("from lada import jobs as ladajobs")
    else: ip.ex("from lada import " + key)


_main()
