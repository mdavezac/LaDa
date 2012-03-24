def listjobs(self, arg):
  """ Lists subjobs. """
  from lada import interactive
  from ..jobs import JobParams
  if interactive.jobdict is None: return
  current = JobParams(jobdict=interactive.jobdict)[interactive.jobdict.name]
  if len(arg) != 0:
    if arg == "all": 
      for job in current.jobdict.root.itervalues():
        if job.is_tagged: continue
        print job.name
      return
    current = current[arg]
  string = ''
  for i, name in enumerate(current.iterkeys()):
    if interactive.jobdict.name == name[:min(len(name), len(interactive.jobdict.name))]:
      string += name[len(interactive.jobdict.name):] + '  '
    else: string += name + '  '
    if (i+1) % 6 == 0: string += '\n'
  print string
