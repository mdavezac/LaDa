def listjobs(self, arg):
  """ Lists sub-folders. """
  from lada import interactive
  from ..jobfolder import JobParams
  if interactive.jobfolder is None: return
  current = JobParams(jobfolder=interactive.jobfolder)[interactive.jobfolder.name]
  if len(arg) != 0:
    if arg == "all": 
      for job in current.jobfolder.root.itervalues():
        if job.is_tagged: continue
        print job.name
      return
    current = current[arg]
  string = ''
  for i, name in enumerate(current.iterkeys()):
    if interactive.jobfolder.name == name[:min(len(name), len(interactive.jobfolder.name))]:
      string += name[len(interactive.jobfolder.name):] + '  '
    else: string += name + '  '
    if (i+1) % 6 == 0: string += '\n'
  print string if i != 0 else "No sub-folders."
