""" Launch interactive function.

 
    This launch strategy will interactively compute each lada job. This will
    block the interpreter.
""" 
__docformat__ = "restructuredtext en"

def launch(self, event, jobfolders):
  """ Launch jobs interactively.

      This call will block until each job is finished in turn.
  """
  from os.path import join, dirname
  from copy import deepcopy
  from ... import default_comm
  
  kwargs = dict(event.kwargs) if event.kwargs is not None else dict()
  comm = deepcopy(default_comm)
  comm['n'] = event.nbprocs
  comm["ppn"] = event.ppn
  kwargs['comm'] = comm


  for current, path in jobfolders:
    # start computations.
    for job in current.itervalues(): 
      name = str(job.name)
      if name[0] == '/': name = name[1:]
      if hasattr(job.functional, 'Extract') and not event.force: 
        p = join(dirname(path), name)
        extract = job.functional.Extract(p)
        if extract.success:
          print "Job {0} completed successfully. It will not be relaunched.".format(name)
          continue
      print "Working on {0} in {1}.".format(name, path)
      kwargs["outdir"] = join(dirname(path), name)
      job.compute(**kwargs)


def completer(self, event, data):
  """ Completer for scattered launcher. """
  from .. import jobfolder_file_completer
  if data[-1] == "--kwargs":
    return [u for u in self.api.user_ns if u[0] != '_' and isinstance(self.api.user_ns[u], dict)]
  elif data[-1] == "--nbprocs": return ['']
  elif data[-1] == "--ppn": return ['']
  result = ['--force', '--kwargs', '--help', '--nbprocs', '--ppn']
  result.extend(jobfolder_file_completer(self, [event.symbol]))
  result = list(set(result) - set(data))
  return result

def parser(self, subparsers, opalls):
  """ Adds subparser for interactive. """ 
  from ... import default_comm
  result = subparsers.add_parser( 'interactive',
                                  description="Launches calculations interactively.\n"\
                                              "Each job will launched one after the other. "\
                                              "This call is *blocking*.",
                                  parents=[opalls] )
  result.add_argument( '--kwargs', type=dict, default={}, dest="kwargs",
                       help="Dictionary which contains arguments for the functionals. "\
                            "\"outdir\" and \"comm\" are added automatically. "\
                            "The functional must accept these arguments." )
  result.add_argument( '--nbprocs', type=int, default=default_comm.get('n', 1),
                       nargs='?', help="Number of processes over which to launch external calculations. "\
                                       "Defaults to {0}.".format(default_comm.get('n', 1)))
  result.add_argument( '--ppn', dest="ppn", default=default_comm.get('ppn', 1), type=int,
                       help="Number of processes per node with which to launch external calculations. "\
                            "Defaults to {0}.".format(default_comm.get('ppn', 1)) )
 
  result.set_defaults(func=launch)
  return result
