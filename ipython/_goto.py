""" IPython goto and iterate magic functions. """
from contextlib  import contextmanager

def goto(self, arg):
  """ Moves current dictionary position and working directory (if appropriate). """
  from os import chdir
  from os.path import exists, join, split as splitpath, isdir
  from . import _get_current_job_params
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
  elif args[0] == "pbs":
   if path == None: 
     ip.user_ns["_lada_error"] = "Cannot go to pbs dir: default dictionary path not set."\
                                 "\nPlease user \"savejobs\"."
     print ip.user_ns["_lada_error"]
     return
   elif not exists(path + ".pbs"):
     ip.user_ns["_lada_error"] = "pbs dir %s does not exist." % (path + ".pbs")
     print ip.user_ns["_lada_error"]
   elif not isdir(path + ".pbs"):
     ip.user_ns["_lada_error"] = "pbs dir %s exists but is not a directory." % (path + ".pbs")
     print ip.user_ns["_lada_error"]
   chdir(path+".pbs")
   return 

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


def iterate(self, event):
  """ Goes to next (untagged) job. """
  from . import _get_current_job_params
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
  from os.path import exists, isdir
  import IPython
  from . import _get_current_job_params
  ip = self.api
  current, path = _get_current_job_params(self, 0)
  if current == None: raise IPython.ipapi.TryNext
  if len(event.line.split()) > 2: raise IPython.ipapi.TryNext

  has_pbs = False
  if path != None:
    if exists(path + ".pbs") and isdir(path + ".pbs"):
      has_pbs = True

  if '/' in event.symbol:
    subkey = ""
    for key in event.symbol.split('/')[:-1]: subkey += key + "/"
    try: subdict = current[subkey]
    except KeyError: raise IPython.ipapi.TryNext
    if hasattr(subdict, "children"): 
      if hasattr(subdict.children, "keys"):
        return [subkey + a + "/" for a in subdict.children.keys()]
    raise IPython.ipapi.TryNext
  else:
    result = [a + "/" for a in current.children.keys()]
    result.extend(["/", "next", "reset"])
    if has_pbs: result.append("pbs")
    if current.parent != None: result.append("../")
    if "_lada_subjob_iterated" in ip.user_ns:
      if len(ip.user_ns["_lada_subjob_iterated"]): result.append("previous")
    return result

