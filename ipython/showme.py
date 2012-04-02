""" IPython showme magic function. """

def showme(self, event):
  """ Edits job parameters, view files. """
  from lada import interactive
  # gets dictionary, path.
  if interactive.jobdict is None:
    print "No current jobdictionary."
    return
  # splits argumetns, removes decorative keywords.
  args = [u for u in event.split() if u not in ["in"]]
  # nothing to do.
  if len(args) == 0: 
    print "Showme requires an argument."
    print "Find out which via tab completion."
    return
  elif len(args) > 1: 
    print "Showme accepts only one argument."
    return

  if args[0] in ["pbserr", "pbsout", "pbs"]:
    showme_pbs(self, args[0])
  elif args[0] == 'functional': 
    showme_functional(self)
  elif args[0] in interactive.jobdict.params:
    showme_param(self, args[0])
  else: print "Unknown job parameter {0}.".format(args[0])

def showme_functional(self):
  """ Edits functional for current job. """
  from types import FunctionType
  from tempfile import NamedTemporaryFile
  from os import remove, stat
  from lada import interactive
  from ..misc import read_input
  if interactive.jobdict.functional is None:
    print "No current functional."
    print "Please first set it with jobparams."
    return
  try: # try/finally section will removed namedtemporaryfile.
    # want .py suffix to get syntax highlighting in editors.
    with NamedTemporaryFile("w", delete = False, suffix='*.py') as file:
      filename = file.name
      # editing INCAR.
      if isinstance(interactive.jobdict.functional, FunctionType):
        file.write('from {0.__module__} import {0.__name__}\n'\
                   'functional = {0.__name__}'\
                   .format(interactive.jobdict.functional))
      else: 
        string = repr(interactive.jobdict.functional)
        if len(string) > 1 and string[0] == '<' and string[-1] == '>':
          print "Functional cannot be represented."
          print "Please use jobparams to modify it."
          return
        file.write(string)

    # lets user edit stuff.
    time0 = stat(filename)[-2]
    self.magic("edit -x {0}".format(filename))
    if stat(filename)[-2] == time0: return

    # change jobparameters.
    input = read_input(filename)
    interactive.jobdict.functional = input.functional
  finally:
    try: remove(filename)
    except: pass

def showme_param(self, arg):
  """ Edits a job parameter. """
  from tempfile import NamedTemporaryFile
  from os import remove, stat
  from re import search, M
  from lada import interactive
  from ..misc import read_input

  if arg not in interactive.jobdict.params: 
    print "{0} is not a jobparameter.".format(arg)
    return 

  try: # try/finally section will removed namedtemporaryfile.
    # want .py suffix to get syntax highlighting in editors.
    with NamedTemporaryFile("w", delete = False, suffix='*.py') as file:
      filename = file.name
      string = repr(interactive.jobdict.params[arg])
      if len(string) > 1 and string[0] == '<' and string[-1] == '>':
        print "Parameter {0} cannot be represented.".format(arg)
        print "Please use jobparams to modify it."
        return
      if interactive.jobdict.params[arg].__class__.__module__ != '__builtin__':
        file.write( 'from {0.__class__.__module__} import {0.__class__.__name__}\n'\
                    .format(interactive.jobdict.params[arg]))
        obre = search( '\s*(\S+)\s*=\s*{0.__class__.__name__}\s*\('\
                       .format(interactive.jobdict.params[arg]),
                       string, M )
      else: obre = None
      if obre is None: 
        othername = arg
        string = string.replace('\n', '\n' + ''.join([' ']*(len(arg)+3)))
        file.write('{0} = {1}'.format(arg, string))
      else: 
        othername = obre.group(1)
        file.write(string)
      
    # lets user edit stuff.
    time0 = stat(filename)[-2]
    self.magic("edit -x {0}".format(filename))
    if stat(filename)[-2] == time0: return
  
    # change jobparameters.
    input = read_input(filename)
    if not hasattr(input, othername): 
      print "Cannot find {0} in file. Aborting.".format(othername)
      return
    interactive.jobdict.params[arg] = getattr(input, othername) 
  finally:
    try: remove(filename)
    except: pass

def showme_pbs(self, which):
  """ Shows pbs files through less. """
  from os.path import join, exists
  from glob  import glob
  from operator import itemgetter
  from .. import interactive

  if not interactive.jobdict.is_job:
    print "Current position in job-dictionary is not a job."
    return
  filename = interactive.jobdict.name[:-1].replace("/", ".")
  if which in ["pbserr", "pbsout"]:
    prefix = "out" if which == "pbsout" else "err"
    filename = join(interactive.jobdict_path + ".pbs", prefix + filename)
    filenames = glob(filename+'.*')
    numbers   = [(i, int(u[len(filename)+1:].split('.')[-1])) for i, u in enumerate(filenames)]
    if len(numbers) == 0: filename = None
    else: filename = filenames[ sorted(numbers, key=itemgetter(1))[-1][0] ]
  else: filename = join(interactive.jobdict_path + ".pbs", filename[1:] + ".pbs")
  
  if filename is None or (not exists(filename)): 
    print "Could not find {0}({1}).".format(which, filename)
    return

  self.system("less {0}".format(filename))

def completer(self, event):
  """ Completer for showme. """
  from os.path import exists, join, dirname
  from glob import glob
  from lada import interactive
  if interactive.jobdict is None: return ['']
  if interactive.jobdict.functional is None: return ['']
  result = ['functional'] + interactive.jobdict.params.keys()
  if interactive.jobdict.is_job and interactive.jobdict_path is not None: 
    jobname = interactive.jobdict.name[:-1].replace("/", ".")
    filename = join(interactive.jobdict_path + ".pbs", 'err' + jobname)
    if len(glob(filename + ".*")) > 0: result.append('pbserr')
    filename = join(interactive.jobdict_path + ".pbs", 'out' + jobname)
    if len(glob(filename + ".*")) > 0: result.append('pbsout')
    filename = join(interactive.jobdict_path + ".pbs", jobname[1:] + ".pbs")
    if exists(filename): result.append('pbs')
  return result
