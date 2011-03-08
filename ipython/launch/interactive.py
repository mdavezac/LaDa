""" Launch interactive function.

 
    This launch strategy will interactively compute each lada job. This will
    block the interpreter.
""" 
__docformat__ = "restructuredtext en"

def launch(self, event, jobdicts):
  """ Launch jobs interactively.

      This call will block until each job is finished in turn.
  """
  import re
  from os.path import join, dirname
  from ... import lada_with_mpi
  from ...opt.changedir import Changedir
  from ..mpi import world
  
  ip = self.api
  ip.user_ns.pop("_lada_error", None)

  kwargs = dict(event.kwargs) if event.kwargs != None else dict()
  kwargs['comm'] = world

  for current, path in jobdicts:
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


def completer(self, info, data, computer):
  """ Completer for scattered launcher. """
  from .._explore import _glob_job_pickles
  if    (len(info.symbol) == 0 and data[-1] == "--kwargs") \
     or (len(info.symbol) > 0  and data[-2] == "--kwargs"):
    return [u for u in iself.api.user_ns if u[0] != '_' and isinstance(self.api.user_ns[u], dict)]
  result = ['--force', '--kwargs', '--help']
  result.extend( _glob_job_pickles(self.api, info.symbol) )
  result = list(set(result) - set(data))
  return result

def parser(self, subparsers, which, opalls):
  """ Adds subparser for interactive. """ 
  result = subparsers.add_parser( 'interactive',
                                  description="Launches calculations interactively.\n"\
                                              "Each job will launched one after the other using "\
                                              "all available processors (unless lada_with_mpi "\
                                              "is false). This call is *blocking*.",
                                  parents=[opalls] )
  result.add_argument( '--kwargs', type=dict, default={}, dest="kwargs",
                       help="Dictionary which contains arguments for the functionals. "\
                            "\"outdir\" and \"comm\" are added automatically. "\
                            "The functional must accept these arguments." )
  result.add_argument( '--force', action="store_true", dest="force", \
                       help="Launches all untagged jobs, even those "\
                            "which completed successfully." )
  result.set_defaults(func=launch)
  return result
