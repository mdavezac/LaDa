import _doc
__doc__ = _doc.__doc__
__docformat__ = "restructuredtext en"

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
    return current, None
  path = ip.user_ns["current_jobdict_path"]
  return current, path

def listjobs(self, arg):
  """ Lists subjobs. """
  ip = self.api
  current, path = _get_current_job_params(self, 1)
  ip.user_ns.pop("_lada_error", None)
  if current is None: return
  if len(arg) != 0:
    if arg == "all": 
      for job in current.root.itervalues():
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
  from ..jobs import JobParams, MassExtract as Collect, save as savejobs
  ip = self.api
  # gets dictionary, path.
  current, path = _get_current_job_params(self, 1)
  ip.user_ns.pop("_lada_error", None)
  if current is None:
    ip.user_ns["_lada_error"] = "No job-dictionary to save."
    print ip.user_ns["_lada_error"] 
    return
  args = [u for u in event.split() ]
  if len(args) == 0: 
    if path is None: 
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
        a = raw_input("File %s already exists.\nOverwrite? [y/n] " % (path))
      if a == 'n':
       ip.user_ns["_lada_error"] = "User said no save."
       return
    savejobs(current.root, path, overwrite=True, timeout=10) 
  elif len(args) == 1:
    if exists(args[0]): 
      if not isfile(args[0]): 
        ip.user_ns["_lada_error"] = "%s is not a file." % (path)
        print ip.user_ns["_lada_error"] 
        return
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File %s already exists.\nOverwrite? [y/n] " % (args[0]))
      if a == 'n':
       ip.user_ns["_lada_error"] = "User said no save."
       return
    savejobs(current.root, args[0], overwrite=True, timeout=10) 
    ip.user_ns["current_jobdict_path"] = abspath(args[0])
    if "collect" not in ip.user_ns: ip.user_ns["collect"] = Collect(dynamic=True)
    if "jobparams" not in ip.user_ns: ip.user_ns["jobparams"] = JobParams()
  else:
    ip.user_ns["_lada_error"] = "Invalid call to saveto."
    print ip.user_ns["_lada_error"] 


def current_jobname(self, arg):
  """ Returns current jobname. """
  ip = self.api
  if "current_jobdict" not in ip.user_ns: return
  print ip.user_ns["current_jobdict"].name
  return

def cancel_completer(self, info):
  return qstat(self, info.symbol).fields(-1)[1:]

def cancel_jobs(self, arg):
  """ Cancel jobs which grep for whatever is in arg.
  
      For instance, the following cancels all jobs with "anti-ferro" in their
      name.

      >>> %cancel_jobs "anti-ferro"
  """
  from lada import lada_with_slurm
  from .qstat import qstat
  arg = str(arg[1:-1])
  if len(arg) != 0: 
    result = qstat(self, arg)
    for u, name in zip(result.fields(0), result.fields(-1)):
      print "cancelling %s." % (name)
    message = "Are you sure you want to cancel the jobs listed above? [y/n] "
  else: message = "Cancel all jobs? [y/n] "
  a = ''
  while a not in ['n', 'y']: a = raw_input(message)
  if a == 'n': return
  
  cmd = "scancel " if lada_with_slurm  else  "qdel "
  result = qstat(self, arg)
  for u, name in zip(result.fields(0), result.fields(-1)): self.api.system(cmd + str(u))


def ipy_init():
  """ Initialises ipython session. 

      In order to initialize an ipython session with the lada interface, the
      following lines can be introduces in the main function of your
      ipy_user_conf.py

      >>> try: import lada.ipython 
      >>> except ImportError:
      >>>   print "Could not find module lada.ipython."
      >>>   pass
      >>> else: lada.ipython.ipy_init()
  """ 
  try: import IPython.ipapi
  except: pass
  else:
    import lada
    from ._goto import goto, goto_completer
    from ._explore import explore, explore_completer
    from ._showme import showme, showme_completer
    from .launch import launch, completer as launch_completer
    from .export import export, completer as export_completer
    from .record import record, completer as record_completer
    from .qstat import qstat
    
    ip = IPython.ipapi.get()
    ip.expose_magic("explore", explore)
    ip.expose_magic("goto", goto)
    ip.expose_magic("listjobs", listjobs)
    ip.expose_magic("showme", showme)
    ip.expose_magic("savejobs", saveto)
    ip.expose_magic("launch", launch)
    ip.expose_magic("qstat", qstat)
    ip.expose_magic("cancel_jobs", cancel_jobs)
    ip.expose_magic("export", export)
    ip.expose_magic("record", record)
    ip.set_hook('complete_command', goto_completer, re_key = '\s*%?goto')
    ip.set_hook('complete_command', showme_completer, re_key = '\s*%?showme')
    ip.set_hook('complete_command', explore_completer, re_key = '\s*%?explore')
    ip.set_hook('complete_command', launch_completer, re_key = '\s*%?launch')
    ip.set_hook('complete_command', cancel_completer, re_key = '\s*%?cancel_jobs')
    ip.set_hook('complete_command', export_completer, re_key = '\s*%?export')
    ip.set_hook('complete_command', record_completer, re_key = '\s*%?record')

    if 'ladabase' in lada.__all__:
      if getattr(lada, 'add_push_magic_function', False): 
        def push(self, cmdl):
          """ Magic function to push to database. 
  
              This detour makes sure ladabase is only loaded on call, 
              thus speeding up load time when starting ipython.
          """
          from lada.ladabase.push import push
          return push(self, cmdl)
        ip.expose_magic("push", push)
      # Don't try and start the pymongo interface unless explicitly requested.
      if lada.ladabase_doconnect: 
        from lada.ladabase import Manager
        try: manager = Manager()
        except: pass
        else: ip.user_ns['ladabase'] = manager

    if len(lada.auto_import_modules) > 0: 
      mods = __import__('lada', globals(), locals(), lada.auto_import_modules)
      for key in set(lada.auto_import_modules) & set(lada.__all__):
        print "Importing", key, "from lada into namespace."
        if key == "jobs": ip.user_ns['ladajobs'] = mods.jobs
        elif key == "ladabase": ip.user_ns['ladabase_module'] = mods.ladabase
        else: ip.user_ns[key] = getattr(mods, key)
    else:
      print "Lada modules were not loaded into user namespace."
      print "Lada interface (explore...) is available."
