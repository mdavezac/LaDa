def listfolders(self, arg):
  """ Lists sub-folders. """
  from fnmatch import fnmatch
  from pylada import interactive
  from ..jobfolder import JobParams
  if interactive.jobfolder is None: return
  if len(arg) == 0:
    string = ''
    for i, name in enumerate(interactive.jobfolder.children.iterkeys()):
      string += name + '  '
      if (i+1) % 6 == 0: string += '\n'
    print string if len(string) != 0 else "No sub-folders."
    return
  elif 'all' in arg.split():
    current = JobParams(jobfolder=interactive.jobfolder)[interactive.jobfolder.name]
    for job in current.jobfolder.root.itervalues():
      if job.is_tagged: continue
      print job.name
    return
  else:
    dirs = arg.split('/')
    result = set()
    for name in interactive.jobfolder.iterleaves():
      name = name[len(interactive.jobfolder.name):]
      if len(name) == 0: continue
      names = name.split('/')
      if len(names) < len(dirs): continue
      if all(fnmatch(u,v) for u,v in zip(names, dirs)): 
        result.add('/'.join(names[:len(dirs)]))
    for i, string in enumerate(result):
      print string,
      if (i+1) % 6 == 0: print '\n'
