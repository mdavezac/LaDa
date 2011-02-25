""" Launch interactive function.

 
    This launch strategy will interactively compute each lada job. This will
    block the interpreter.
""""
__docformat__ = "restructuredtext en"

def interactive(self, event):
  """ Launch scattered jobs: one job = one pbs script. """
  import re
  from os.path import join
  from ..opt.changedir import Changedir
  from .. import lada_with_mpi
  
  ip = self.api
  ip.user_ns.pop("_lada_error", None)

  # creates list of dictionaries.
  pickles = set(event.pickle) - set([""])
  if len(pickles) > 0: 
    jobdicts = []
    for p in pickles:
      try: d = jobs.load(path=p)
      except: 
        print "JobDict could not be loaded form {0}.".format(p)
        return
      jobdicts.append((d, p))
  else: # current job dictionary.
    current, path = _get_current_job_params(self, 2)
    if current == None: return
    if path == None: return
    # saving pickle
    saveto(self, path)
    if "_lada_error" in ip.user_ns:
      if ip.user_ns["_lada_error"] == "User said no save.":
        print "Job-dictionary not saved = jobs not launched."
      return
    jobdicts = [(current, path)]

  kwargs = dict(event.kwargs) if event.kwargs != None else dict()
  if lada_with_mpi: 
    from boost.mpi import world
    kwargs["comm"] = world
  else: kwargs["comm"] = None

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
      kwargs["outdir"] = join(path, name)
      job.compute(**kwargs)


def completer(self, info, data, computer):
  """ Completer for scattered launcher. """
  from .._explore import _glob_job_pickles
  if    (len(info.symbol) == 0 and data[-1] == "--kwargs") \
     or (len(info.symbol) > 0  and data[-2] == "--kwargs"):
    return [u for u in ip.user_ns if u[0] != '_' and isinstance(ip.user_ns[u], dict)]
  result = ['--kwargs', '--help']
  result = list(set(result) - set(data))
  result.extend( _glob_job_pickles(ip, info.symbol) )
  return result

def parser(self, subparsers, which, opalls):
  """ Adds subparser for interactive. """ 
  result = subparsers.add_parser( 'interactive',
                                  description='Launches calculations interactively.',\
                                  parents=[opalls] )
  result.add_argument( '--kwargs', type=dict, default={}, dest="kwargs",
                       help="Dictionary which contains arguments for the functionals."\
                            "\"outdir\" and \"comm\" are added automatically. "\
                            "The functional must accept these arguments." )
  result.add_argument( '--force', action="store_true", dest="force", \
                       help="Launches all untagged jobs, even those "\
                            "which completed successfully." )
  result.set_defaults(func=launch)
  return result
