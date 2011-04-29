""" Allows to rapidly store and reload information to disk. 

    After a lengthy calculation it may be practical to save the results to file. 

    >>> result = 5 
    >>> record result

    The above save the variable ``result`` to disk.
    It can then be retrieved, at the expense of overwriting the value in the
    user namespace.

    >>> result
    5
    >>> result = 8
    >>> record --load result
    Reloaded result(=5) from record .lada_record.
    >>> result
    5

    All records can be viewed and listed using:

    >>> record --view
    result = 5
    >>> record --list
    ['result']
"""
def record(self, cmdl):
  """ Records variables to file. """
  from pickle import load, dump, dumps
  from argparse import ArgumentParser
  from os.path import exists
  from lada.opt import RelativeDirectory
  
  parser = ArgumentParser( prog="%record",
                           description="Allows to rapidly store and reload information to disk.")
  parser.add_argument( '--file', dest="filename", default=".lada_record", type=str, 
                       help="Name of file where variables are recorded." )
  parser.add_argument( 'vars', metavar='VAR', type=str, nargs='*',
                       help='Name of the variable(s) to record or remove.' )
  group = parser.add_mutually_exclusive_group()
  group.add_argument( '--remove', action="store_true", dest='remove',
                       help="Removes variable from record if present." )
  group.add_argument( '--update', action="store_true", dest='update',
                       help="Adds/updates variables to record. Default." )
  group.add_argument( '--list', action="store_true", dest='list',
                       help="Lists objects in record." )
  group.add_argument( '--view', action="store_true", dest='view',
                       help="View a given record, if possible printable."
                            "If not VAR argument is given, all records are listed." )
  group.add_argument( '--load', action="store_true", dest='load',
                       help="Reloads variables from record into user namespace."\
                            "If no VAR are given, reload all records." )

  try: args = parser.parse_args(cmdl.split())
  except SystemExit:
    if '-h' in cmdl: print __doc__[__doc__.find('\n'):].replace('\n    ', '\n')
    return None
  if len(args.vars) == 0 and not (args.list or args.load or args.view):
    parser.print_help()
    print "\n*****************\n"\
          "At least on VAR argument required.\n"\
          "*****************"
    return 

  # open record file.
  path = RelativeDirectory(args.filename).path
  if exists(path): 
    with open(path, 'r') as file: dummy, values = load(file)
  elif args.remove or args.load or args.list:
    print "Path {0} does not exist.\nNo records yet.\n".format(args.filename)
    return
  else: values = {}

  has_changed = False
  if args.list: return values.keys()
  if args.view: 
    if len(args.vars) == 0: args.vars = values.iterkeys()
    for var in args.vars:
      if var not in values:
        print "Could not find {0} in {1}.".format(var, args.filename)
        continue
      try: string = str(values[var])
      except: 
        print "{0} in record {1} is not printable.".format(var, args.filename)
        continue
      else:
        if len(string) > 30: string = string[:25] + "..."
        if "\n" in string: string = string[:string.find('\n')] + "..."
        print "{0} = {1}".format(var, string)
  elif args.remove: 
    # Remove argumenst from record.
    for var in set(args.vars):
      if var not in values: 
        print "{0} could not be found in record file {1}.".format(var, args.filename)
        continue
      try: string = str(values.pop(var))
      except: print "Removing {0} from record {1}.".format(var, args.filename)
      else:
        if len(string) > 30: string = string[:25] + "..."
        if "\n" in string: string = string[:string.find('\n')] + "..."
        print "Removing {0}(={1}) from record {2}".format(var, string, args.filename)
      has_changed = True
  elif args.load: 
    if len(args.vars) == 0: args.vars = values.iterkeys()
    for key in args.vars:
      if key not in values:
        print "Could not find {0} in {1}.".format(key, args.filename)
        continue
      self.api.user_ns[key] = values[key]
      try:string = str(values[key]) 
      except: print "Reloaded {0} from record {1}.".format(key, args.filename)
      else:
        if len(string) > 30: string = string[25:] + "..."
        if "\n" in string: string = string[:string.find('\n')] + "..."
        print "Reloaded {0}(={2}) from record {1}.".format(key, args.filename, string)
    return
  else:
    # Add arguments to record.
    for key in set(args.vars): 
      # checks value exists in user namespace.
      if key not in self.api.user_ns: 
        print "Could not find {0} in user namespace.".format(key)
        continue
      # checks the argument can be pickled.
      try: dumps(self.api.user_ns[key])
      except:  print "{0} is not pickleable.".format(key)
      else:
        values[key] = self.api.user_ns[key]
        has_changed = True

  if has_changed:
    with open(path, 'w') as file: dump( ("This is a record.", values), file)

         
def completer(self, event): 
  """ Completer for %record magic function. """ 
  from os.path import isdir

  result = []
  data = event.line.split()[1:]
  if '--file' not in data: result.append('--file') 
  if len(set(['--list', '--view', '--load', '--remove', '--update']).intersection(set(data))) == 0:
    result.extend(['--list', '--view', '--load', '--remove'])
  if    (len(event.symbol) == 0 and len(data) > 0 and data[-1] == "--file") \
     or (len(event.symbol) > 0  and len(data) > 1 and data[-2] == "--file"):
    other = event.symbol
    string = '%mglob "cont:This is a record." {0}*'.format(other)
    result = [u for u in self.api.magic(string)]
    string = '%mglob dir:{0}*'.format(other)
    result.extend([u for u in self.api.magic(string)])
    if isdir(other) and other[-1] != '/':
      string = '%mglob "cont:This is a record." {0}/*'.format(other)
      result.extend([u for u in self.api.magic(string)])
  else:
    result.extend([u for u in self.api.user_ns.iterkeys() if u[0] != '_'])
  return result

