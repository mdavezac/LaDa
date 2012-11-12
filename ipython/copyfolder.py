""" Magic function to copy a functional or subtree. """
def copy_folder(self, event):
  """ Copies a jobfolder somewhere. 


      Emulates bash's ``cp`` command. By default, it only copies one job.
      However, it can be made to copy a full tree as well, using the ``-r``
      directive.  When it finds it would overwrite a non-empty jobfolder, it
      prompts the user interactively, unless the '-f' directive is given.


      >>> copyfolder thisjob thatjob

      Copies ``thisjob`` to ``thatjob``, but not the subfolders of ``thisjob``.

      By default, the jobfolders are deepcopied. This means that the none of
      the arguments or functionals of the current job are shared by the copied
      jobs. However, the relation-ships between the current jobs are retained
      in the destinations. In other words, if jobs 'JobA' and 'JobA/JobB' share
      the same 'structure' variable object and they are copied to 'JobC' (with
      '-r'), then two new subfolders are created, 'JobC/JobA'
      and'JobC/JobA/JobB'. Both of these subfolders will share a reference to
      the same `structure` object.  However their `structure` object is
      different from that of the original `JobA` and `JobB`. This feature can
      be turned off with the option `--nodeepcopy` (in which case, `structure`
      would be shared by all source and destination folders). Furthermore, if
      '--nodeepcopy' is used, then the functionals are shared between source
      and destination.
  """
  from argparse import ArgumentParser
  from os.path import join, normpath, relpath
  from copy import deepcopy
  from ..interactive import jobfolder as cjf

  parser = ArgumentParser(prog='%copyfolder',
               description='Copies a jobfolder from one location to another')
  parser.add_argument('-f', '--force', action='store_true', dest='force',
               help='Does not prompt when overwriting a job-folder.')
  parser.add_argument('-r', '--recursive', action='store_true',
               dest='recursive', help='Whether to copy subfolders as well.')
  parser.add_argument('--nodeepcopy', action='store_true',
               help='Destination folders will share the parameters '           \
                    'from the original folders.')
  parser.add_argument('source', type=str, metavar='SOURCE',
               help='Jobfolder to copy')
  parser.add_argument('destination', type=str, metavar='DESTINATION',
               help='Destination folder')

  try: args = parser.parse_args(event.split())
  except SystemExit: return None

  shell = get_ipython()

  # gets current folder.
  if 'jobparams' not in shell.user_ns:
    print 'No jobfolder currently loaded.'
    return 
  jobparams = shell.user_ns['jobparams']

  # normalize destination.
  if args.destination[0] != ['/']:
    destination = normpath(join(cjf.name, args.destination))
    if destination[0] != '/':
      print 'Incorrect destination', destination
      return 
  else: destination = normpath(args.destination)

  # create list of source directories
  if args.source[0] != ['/']:
    source = normpath(join(cjf.name, args.source))
    if source[0] != '/':
      print 'Incorrect source', source
      return 
  else: source = normpath(args.source)
  if source not in cjf: 
    print 'Source', source, 'does not exist'
    return
  if destination == source: 
    print 'Source and destination are the same'
    return 
  rootsource = source
  pairs = []
  if cjf[rootsource].is_executable and not cjf[rootsource].is_tagged:
    pairs = [(source, destination)]
  if args.recursive:
    for source in jobparams[rootsource]:
      if not cjf[source].is_executable: continue
      if cjf[source].is_tagged: continue
      pairs.append((source, join(destination, relpath(source, rootsource))))
  if len(pairs) == 0:
    print "Nothing to copy."
    return 

  # now performs actual copy
  root = deepcopy(cjf.root) if not args.nodeepcopy else cjf.root
  for source, destination in pairs:
    # gets the jobfolder source.
    jobsource = root[source]
    # gets the jobfolder destination.
    jobdest = cjf
    for name in destination.split(): jobdest = jobdest / name
    # something already exists here
    if jobdest.is_executable and not args.force: 
      print 'Copying', jobsource.name, 'to', jobdest.name
      a = ''
      while a not in ['n', 'y']:
        a = raw_input( '{0} already exists. Overwrite? [y/n]'
                       .format(jobdest.name) )
      if a == 'n': 
        print jobdest.name, 'not overwritten.'
        continue
    # now copies folder items.
    for key, value in jobdest.__dict__.iteritems(): 
      if key not in ['children', 'parent', 'param']:
        jobdest.__dict__[key] = value
    jobdest.params = jobsource.params.copy()

def completer(self, event):
  from IPython import TryNext
  from lada import interactive

  line = event.line.split()[1:]
  result = []
  if '-h' not in line and '--help' not in line: result.append('-h')
  if '-r' not in line and '--recursive' not in line: result.append('-r')

  noop = [k for k in line if k[0] != '-']
  if len(noop) > 2: raise TryNext()

  if '/' in event.symbol:
    subkey = ""
    for key in event.symbol.split('/')[:-1]: subkey += key + "/"
    try: subdict = interactive.jobfolder[subkey]
    except KeyError: raise TryNext
    if hasattr(subdict, "children"): 
      if hasattr(subdict.children, "keys"):
        return [subkey + a + "/" for a in subdict.children.keys()]
    raise TryNext()
  else:
    result += [a + "/" for a in interactive.jobfolder.children.keys()]
    result.extend(["/"])
    if interactive.jobfolder.parent is not None: result.append("../")
    if len(getattr(interactive, "_lada_subjob_iterated", [])) != 0:
      result.append("previous")
    return result
