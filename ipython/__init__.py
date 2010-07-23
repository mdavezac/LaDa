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
  from ._launch import launch, launch_completer

  ip = IPython.ipapi.get()
  ip.expose_magic("explore", explore)
  ip.expose_magic("goto", goto)
  ip.expose_magic("listjobs", listjobs)
  ip.expose_magic("jobname", current_jobname)
  ip.expose_magic("iterate", iterate)
  ip.expose_magic("showme", showme)
  ip.expose_magic("savejobs", saveto)
  ip.expose_magic("fakerun", fakerun)
  ip.expose_magic("launch_scattered", launch_scattered)
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
