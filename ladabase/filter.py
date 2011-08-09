""" IPython interface to filter through the database. """

def filters(self, cmdl):
  """ Provides a stack of filters through which to filter the database. """
  import argparse
  import lada
  from .massextract import MassExtract

  if not hasattr(lada, "filters_stack"):
    lada.filters_stack = [("True", MassExtract(namespace=ip.user_ns))]
  stack = lada.filters_stack

  # empty commandline. Print filters.
  if len(cmdl.replace(' ', '')): cmld = "--view"

  parser = argparse.ArgumentParser(prog='%filter',
                     formatter_class=argparse.RawDescriptionHelpFormatter,
                     epilog = \
"""
Creates stack of filters to navigate database. Only those elements satisfying
all filters can be analyzed with "collect". Each call to this magic function
will add a further filter. Filters are specified as below:

  >>> %filters 'Ba' in extract.species

Where *extract* is an object used to access job-parameters from the outcar.
To view which parameters are accessible, type:

  >>> collect.[TAB]

including the \".\" (dot) and hit TAB.
""" )

  group = parser.add_mutually_exclusive_group()
  group.add_argument( 'exprs', type=str, nargs="*", default="",
                       help='Python expression returning True or False.')
  group.add_argument( '--iterate', action="store_true", 
                       dest="iterate", help='Iterates over current job.' )
  group.add_argument( '--view', action="store_true", 
                       dest="view", help='Prints all filters.' )
  group.add_argument( '--pop', action="store_true", 
                       dest="pop", help='Remove last filter.' )
  group.add_argument( '--outcar', action="store_true", 
                       dest="outcar", help='View current OUTCAR, if there is only one.' )

  try: args = parser.parse_args(cmdl.split())
  except SystemExit as e:
    
    return None

  if args.view: 
    if len(stack) == 1:
      print "No filters yet."
    else:
      print "Filters:"
      for i, (filter, value) in enumerate(stack):
        if i == 0: continue
        else: print "{0:>2}. {1}".format(i, filter)
  elif args.iterate: 
    if not hasattr(lada, "iterstack"): lada.iterstack = 0
    else: stack.pop(-1)
    stopiter = True
    for i, key in enumerate(stack[-1][1]):
      if i == lada.iterstack: stopiter = False; break
    if stopiter == True:
      print "Reached last filtered element."
      key = stack[-1][1].keys()[-1]
    else: lada.iterstack += 1
    stack.append(("id == {0}".format(repr(key)), stack[-1][1]["id == {0}".format(repr(key))]))
  elif args.pop and len(stack) > 1:
    if hasattr(lada, "iterstack"): del lada.iterstack
    stack.pop(-1)
  elif args.outcar:
    if len(stack[-1][1]) != 1:
      print "outcar only valid if database is filtered to a single outcar."
    else: 
      from tempfile import NamedTemporaryFile
      try: 
        filename = None
        with NamedTemporaryFile("w", delete=False) as file:
          filename = file.name
          with stack[-1][1].values()[0].__outcar__() as foutcar:
            file.write(foutcar.read())
        self.api.system("less {0}\n".format(filename))
      except: raise
      finally:
        from os import remove
        if filename != None:
         try: remove(filename)
         except: pass
  else: 
    if hasattr(lada, "iterstack"): del lada.iterstack
    cmdl = " ".join(args.exprs)
    try:
      add = stack[-1][1][cmdl]
      print "Found {0} calculations statisfying the filters.".format(len(add))
    except: 
      print "Could not make sense of expression {0}.".format(cmdl)
      if len(stack[-1][1]): 
        self.api.user_ns["extract"] = stack[-1][1].values()[0]
        print "I have put an \"extract\" variable in your namespace for you to try it out."
      else: print "But in any case, there is no jobs to filter anymore."
    else:
      if len(add) == 0: 
        print "This filter has not been added to the list."
      else: stack.append((cmdl, add))

  lada.filters_stack = stack
  self.api.user_ns["collect"] = stack[-1][1]

def completer(self, event):
  """ Completer for filter magic function. """
  return ["--iterate", "--view", "--pop"]

def init(ip):
  """ Initializes filters stuff. """
  from .massextract import MassExtract
  import lada
  ip.expose_magic("filters", filters)
  ip.set_hook('complete_command', completer, re_key = '\s*%?filter')
  if not hasattr(lada, "filters_stack"):
    lada.filters_stack = [("True", MassExtract(namespace=ip.user_ns))]
  ip.user_ns["collect"] = lada.filters_stack[-1][1]

