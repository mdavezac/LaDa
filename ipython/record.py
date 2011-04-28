""" Allows to rapidly store information to disk. """
def record(self, cmdl):
  """ Records variables to file. """
  from pickle import load, dump, dumps
  from argparse import ArgumentParser
  from os.path import exists
  from lada.opt import RelativeDirectory
  
  parser = ArgumentParser( prog="%record", 
                     description="Adds, removes, or updates variable to record file." )
  parser.add_argument( '--file', dest="filename", default=".lada_record", type=str, nargs='?', 
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
  group.add_argument( '--load', action="store_true", dest='load',
                       help="Reloads variables from record into user namespace." )


  try: args = parser.parse_args(cmdl.split())
  except SystemExit: return None

  # open record file.
  path = RelativeDirectory(args.filename).path
  if exists(path): 
    with open(path, 'r') as file: dummy, values = load(file)
  elif args.remove:
    print "Path {0} does not exist.\n".format(args.filename)
    return

  has_changed, result = False, None
  if args.list: result = values.keys()
  if args.remove: 
    # Remove argumenst from record.
    for var in set(args.vars):
      if var not in values: 
        print "{0} could not be found in record file {1}.".format(var, args.filename)
        continue
      values.pop(var)
      has_changed = True
  elif args.load: 
    self.api.update(values)
    for key in values.iterkeys():
      print "Reloaded {0} from record {1}.".format(key, args.filename)
    return
  else:
    # Add arguments to record.
    for key in set(args.vars): 
      # checks value exists in user namespace.
      if key not in self.api.user_ns: 
        print "Could not find {0} in user namespace.".format(key)
        continue
      # checks the argument can be pickled.
      try: dumps(ip.user_ns[key])
      except:  print "{0} is not pickleable.".format(key)
      else:
        values[key] = self.api.user_ns[key]
        has_changed = True

  if has_changed:
    with open(path, 'w') as file: dump( ("This is a record.", values), file)

  return result

       
         
def completer(self, event): 
  """ Completer for %record magic function. """ 
  from os.path import isdir

  result = []
  data = event.line.split()[1:]
  if '--file' not in data: result.append('--file') 
  if len(set(['--list', '--load', '--remove', '--update']).intersect(set(data))) == 0:
    result.extend(['--list', '--load', '--remove'])
  if len(data) >0 and data[-1] == '--list':
    other = event.symbol
    string = '%mglob "cont:This is a record." {0}*'.format(other)
    result = [u for u in self.api.magic(string)]
    string = '%mglob dir:{0}*'.format(other)
    result.extend([u for u in self.api.magic(string)])
    if isdir(other) and other[-1] != '/':
      string = '%mglob "cont:This is a record." {0}/*'.format(other)
      result.extend([u for u in self.api.magic(string)])
  else: result.extend([u for u in self.api.user_ns.iterkeys() if u[0] != '_'])
  return result
