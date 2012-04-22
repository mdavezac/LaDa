""" IPython showme magic function. """

def showme(self, event):
  """ Edits job parameters, view files. """
  from lada import interactive
  # gets dictionary, path.
  if interactive.jobfolder is None:
    print "No current job-folder."
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

  if args[0] == "params": showme_params(self)
  elif args[0] in ["pbserr", "pbsout", "pbs"]:
    showme_pbs(self, args[0])
  elif args[0] == 'functional': 
    showme_functional(self)
  elif args[0] in interactive.jobfolder.params:
    showme_param(self, args[0])
  else: print "Unknown job parameter {0}.".format(args[0])

def showme_functional(self):
  """ Edits functional for current job. """
  from types import FunctionType
  from tempfile import NamedTemporaryFile
  from os import remove, stat
  from lada import interactive
  from ..misc import read_input
  if interactive.jobfolder.functional is None:
    print "No current functional."
    print "Please first set it with jobparams."
    return
  try: # try/finally section will removed namedtemporaryfile.
    # want .py suffix to get syntax highlighting in editors.
    with NamedTemporaryFile("w", delete = False, suffix='*.py') as file:
      filename = file.name
      # editing INCAR.
      if isinstance(interactive.jobfolder.functional, FunctionType):
        file.write('from {0.__module__} import {0.__name__}\n'\
                   'functional = {0.__name__}'\
                   .format(interactive.jobfolder.functional))
      else: 
        string = repr(interactive.jobfolder.functional)
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
    interactive.jobfolder.functional = input.functional
  finally:
    try: remove(filename)
    except: pass

def showme_param(self, arg):
  """ Edits a job parameter. """
  from tempfile import NamedTemporaryFile
  from os import remove, stat
  from re import search, M, sub
  from lada import interactive
  from ..misc import read_input, import_dictionary, import_header_string

  if arg not in interactive.jobfolder.params: 
    print "{0} is not a jobparameter.".format(arg)
    return 

  try: # try/finally section will removed namedtemporaryfile.
    # want .py suffix to get syntax highlighting in editors.
    with NamedTemporaryFile("w", delete = False, suffix='*.py') as file:
      filename = file.name
      string = repr(interactive.jobfolder.params[arg])
      if len(string) > 1 and string[0] == '<' and string[-1] == '>':
        print "Parameter {0} cannot be represented.".format(arg)
        print "Please use jobparams to modify it."
        return
      if interactive.jobfolder.params[arg].__class__.__module__ != '__builtin__':
        obre = search( r'\s*(\S+)\s*=\s*{0.__class__.__name__}\s*\('\
                       .format(interactive.jobfolder.params[arg]),
                       string, M )
      else: obre = None
      mods = import_dictionary(interactive.jobfolder.params[arg])
      if len(mods) != 0: file.write(import_header_string(mods))
      if obre is None: 
        othername = arg
        string = string.replace('\n', '\n' + ''.join([' ']*(len(arg)+3)))
        file.write('{0} = {1}'.format(arg, string))
      else: 
        string = sub(r"\b{0}\b".format(obre.group(1)), arg, string)
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
    interactive.jobfolder.params[arg] = getattr(input, othername) 
  finally:
    try: remove(filename)
    except: pass

def showme_params(self):
  """ Edits a job parameter. """
  from tempfile import NamedTemporaryFile
  from os import remove, stat
  from re import search, M
  from lada import interactive
  from ..misc import read_input, import_dictionary, import_header_string

  try: # try/finally section will removed namedtemporaryfile.
    # want .py suffix to get syntax highlighting in editors.
    with NamedTemporaryFile("w", delete = False, suffix='*.py') as file:
      filename = file.name
      # get header info.
      mods = {}
      for arg, value in interactive.jobfolder.params.iteritems():
        mods = import_dictionary(value, mods)
        string = repr(value)
        if len(string) > 1 and string[0] == '<' and string[-1] == '>':
          print "Parameter {0} cannot be represented.".format(arg)
          print "Please use jobparams to modify it."
          return
        import_dictionary(value, mods)

      print mods
      if len(mods) != 0: file.write(import_header_string(mods))
      file.write('params = {}\n')

      # now write each parameter in turn
      for arg, value in interactive.jobfolder.params.iteritems():
        string = repr(value)
        if interactive.jobfolder.params[arg].__class__.__module__ != '__builtin__':
          obre = search( '\s*(\S+)\s*=\s*{0.__class__.__name__}\s*\('\
                         .format(interactive.jobfolder.params[arg]),
                         string, M )
        else: obre = None
        if obre is None: 
          othername = arg
          string = string.replace('\n', '\n' + ''.join([' ']*(len(arg)+13)))
          file.write("params['{0}'] = {1}\n".format(arg, string))
        else: 
          file.write(string)
          file.write("params['{0}'] = {1}\n".format(arg, obre.group(1)))
      
    # lets user edit stuff.
    time0 = stat(filename)[-2]
    self.magic("edit -x {0}".format(filename))
    if stat(filename)[-2] == time0: return
  
    # change jobparameters.
    input = read_input(filename)
    if not hasattr(input, 'params'): 
      print "Cannot find {0} in file. Aborting.".format('params')
      return
    interactive.jobfolder.params = getattr(input, 'params') 
  finally:
    try: remove(filename)
    except: pass

def showme_pbs(self, which):
  """ Shows pbs files through less. """
  from os.path import join, exists
  from glob  import glob
  from operator import itemgetter
  from .. import interactive

  if not interactive.jobfolder.is_job:
    print "Current position in job-dictionary is not a job."
    return
  filename = interactive.jobfolder.name[:-1].replace("/", ".")
  if which in ["pbserr", "pbsout"]:
    prefix = "out" if which == "pbsout" else "err"
    filename = join(interactive.jobfolder_path + ".pbs", prefix + filename)
    filenames = glob(filename+'.*')
    numbers   = [(i, int(u[len(filename)+1:].split('.')[-1])) for i, u in enumerate(filenames)]
    if len(numbers) == 0: filename = None
    else: filename = filenames[ sorted(numbers, key=itemgetter(1))[-1][0] ]
  else: filename = join(interactive.jobfolder_path + ".pbs", filename[1:] + ".pbs")
  
  if filename is None or (not exists(filename)): 
    print "Could not find {0}({1}).".format(which, filename)
    return

  self.system("less {0}".format(filename))

def completer(self, event):
  """ Completer for showme. """
  from os.path import exists, join
  from glob import glob
  from lada import interactive
  if interactive.jobfolder is None: return ['']
  if not interactive.jobfolder.is_executable: return ['']
  result = ['functional', 'params'] + interactive.jobfolder.params.keys()
  if interactive.jobfolder_path is not None: 
    jobname = interactive.jobfolder.name[:-1].replace("/", ".")
    filename = join(interactive.jobfolder_path + ".pbs", 'err' + jobname)
    if len(glob(filename + ".*")) > 0: result.append('pbserr')
    filename = join(interactive.jobfolder_path + ".pbs", 'out' + jobname)
    if len(glob(filename + ".*")) > 0: result.append('pbsout')
    filename = join(interactive.jobfolder_path + ".pbs", jobname[1:] + ".pbs")
    if exists(filename): result.append('pbs')
  return result
