""" IPython goto and iterate magic functions. """
from contextlib  import contextmanager

def goto(self, cmdl):
  """ Moves current dictionary position and working directory (if appropriate). """
  from os import chdir
  from os.path import exists, join, split as splitpath, isdir
  from lada import interactive

  if interactive.jobdict is None: 
    print "No current jobs."
    return
  if len(cmdl.split()) == 0:
    if interactive.jobdict_path is None:
      print "Current position in job dictionary:", interactive.jobdict.name
    else:
      print "Current position in job dictionary:", interactive.jobdict.name
      print "Filename of jobdictionary: ", interactive.jobdict_path
    return
  args = cmdl.split()
  if len(args) > 1:
    ip.user_ns["_lada_error"] = "Invalid argument to goto %s." % (cmdl)
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
    if interactive.jobdict_path is None: 
      print "Cannot go to pbs dir: default dictionary path not set."\
            "\nPlease user \"savejobs\"."
      return
    elif not exists(interactive.jobdict_path + ".pbs"):
      print "pbs dir {0} does not exist.".format(interactive.jobdict_path + ".pbs")
      return
    elif not isdir(interactive.jobdict_path + ".pbs"):
      print "pbs dir {0} exists but is not a directory.".format(interactive.jobdict_path + ".pbs")
      return
    chdir(interactive.jobdict_path+".pbs")
    return 

  # case for which precise location is given.
  try: result = interactive.jobdict[args[0]] 
  except KeyError as e: 
    print e
    return 

  interactive.jobdict = result
  good = not interactive.jobdict.is_tagged
  if good:
    for value in interactive.jobdict.values():
      good = not value.is_tagged
      if good: break
  if not good:
    print '**** Current job and subjobs are all off; '\
          'jobparams (except onoff) and collect will not work.'
    return
  if interactive.jobdict_path is None: return
  dir = join(splitpath(interactive.jobdict_path)[0], interactive.jobdict.name[1:]) 
  if exists(dir): chdir(dir)
  return


def iterate(self, event):
  """ Goes to next (untagged) job. """
  from lada import interactive
  if interactive.jobdict is None: return

  args = event.split()
  if len(args) > 1: 
    print "Invalid argument {0}.".format(event)
    return
  elif len(args) == 0:
    if "_lada_subjob_iterator" in interactive.__dict__:
      iterator = interactive._lada_subjob_iterator
    else:
      iterator = interactive.jobdict.root.itervalues()
    while True:
      try: job = iterator.next()
      except StopIteration: 
        print "Reached end of job list."
        return 
      if job.is_tagged: continue
      break
    interactive._lada_subjob_iterator = iterator
    if "_lada_subjob_iterated" not in interactive.__dict__:
      interactive._lada_subjob_iterated = []
    interactive._lada_subjob_iterated.append(interactive.jobdict.name)
    goto(self, job.name)
    print "In job ", interactive.jobdict.name
  elif args[0] == "reset" or args[0] == "restart":
    # remove any prior iterator stuff.
    interactive.__dict__.pop("_lada_subjob_iterator", None)
    interactive.__dict__.pop("_lada_subjob_iterated", None)
    print "In job ", interactive.jobdict.name
  elif args[0] == "back" or args[0] == "previous":
    if "_lada_subjob_iterated" not in interactive.__dict__:
      print "No previous job to go to. "
    else:
      goto(self, interactive._lada_subjob_iterated.pop(-1))
      if len(interactive._lada_subjob_iterated) == 0:
        del interactive._lada_subjob_iterated
      print "In job ", interactive.jobdict.name


def completer(self, event):
  from os.path import exists, isdir
  from IPython import TryNext
  from lada import interactive
  if len(event.line.split()) > 2: raise TryNext

  has_pbs =      interactive.jobdict_path is not None \
             and exists("{0}.pbs".format(interactive.jobdict_path)) \
             and isdir("{0}.pbs".format(interactive.jobdict_path)) 

  if '/' in event.symbol:
    subkey = ""
    for key in event.symbol.split('/')[:-1]: subkey += key + "/"
    try: subdict = interactive.jobdict[subkey]
    except KeyError: raise TryNext
    if hasattr(subdict, "children"): 
      if hasattr(subdict.children, "keys"):
        return [subkey + a + "/" for a in subdict.children.keys()]
    raise TryNext
  else:
    result = [a + "/" for a in interactive.jobdict.children.keys()]
    result.extend(["/", "next", "reset"])
    if has_pbs: result.append("pbs")
    if interactive.jobdict.parent is not None: result.append("../")
    if len(getattr(interactive, "_lada_subjob_iterated", [])) != 0:
      result.append("previous")
    return result

