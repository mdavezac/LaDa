""" IPython showme magic function. """

def showme(self, event):
  """ Edits functional and/or structure. """
  from os import remove, stat
  from os.path import join
  from tempfile import NamedTemporaryFile
  from ..opt import read_input
  from ..crystal import write_poscar, read_poscar, Structure
  from . import _get_current_job_params
  from . import represent_structure_with_POSCAR
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
    if args[0] in ["pbserr", "pbsout", "pbs"]: _showme_pbs(self, args[0]); return
    filename = None
    try: # try/finally section will removed namedtemporaryfile.
      # want .py suffix to get syntax highlighting in editors.
      suffix = ".py" if arg.lower() == "functional" else ''
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
              if key == "functional": continue
              if key == "structure": continue
              file.write("jobparams[\"%s\"] = %s\n" % (key, repr(value)))
        # editing POSCAR.
        elif arg.lower() == "structure":
          # always as vasp5. Makes it easier.
          structure = Structure() if "structure" not in current.jobparams \
                      else current.jobparams["structure"]
          if represent_structure_with_POSCAR:
            write_poscar(current.structure, file, True)
          else: file.write(repr(structure))
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
        if "structure" in current.jobparams: input.jobparams["structure"] = current.structure
        current.functional = input.functional
        current.jobparams = input.jobparams
      elif arg.lower() == "structure": 
        if represent_structure_with_POSCAR:
          d = {}
          if    "structure" in current.jobparams \
             and hasattr(current.jobparams["structure"], __dict__):
            d = current.jobparams["structure"].__dict__
          current.jobparams["structure"] = read_poscar(path=filename)
          current.jobparams["structure"].__dict__.update(d)
        else: 
          globals = ip.user_ns.copy()
          try: execfile(filename, globals)
          except:
            print "Error when executing structure script."
            raise
          else:
            if "structure" not in globals:
              print "Could not find structure in structure script."
              return
            current.jobparams["structure"] = globals["structure"]
            
    finally:
      if filename != None:
        try: remove(filename)
        except: pass

def _showme_pbs(self, which):
  from os.path import join, exists
  from glob  import glob
  from operator import itemgetter
  from . import _get_current_job_params

  ip = self.api
  current, path = _get_current_job_params(self, 1)
  filename = current.name.replace("/", ".")
  if which in ["pbserr", "pbsout"]:
    prefix = "out" if which == "pbsout" else "err"
    filename = join(path + ".pbs", prefix + filename)
    filenames = [u for u in glob(filename+'.*')]
    numbers   = [(i, int(u[len(filename)+1:].split('.')[0])) for i, u in enumerate(filenames)]
    if len(numbers) == 0: filename = None
    else: filename = filenames[ sorted(numbers, key=itemgetter(1))[-1][0] ]
  else: filename = join(path + ".pbs", filename[1:] + ".pbs")
  
  if filename == None or (not exists(filename)):  print "Could not find {0}.".format(which); return

  ip.system("less {0}".format(filename))

def showme_completer(self, event):
  import IPython
  ip = self.api
  if "current_jobdict" not in ip.user_ns: raise IPython.ipapi.TryNext
  return ["structure", "functional", "pbserr", "pbsout", "pbs"]

