""" IPython launch magic function. """

__docformat__ = "restructuredtext en"

def launch_scattered(self, event):
  """ Launch scattered jobs: one job = one pbs script. """
  import re
  from os import environ
  from os.path import split as splitpath, join, exists, abspath, dirname
  from ..opt.changedir import Changedir
  from ..jobs.templates import default_pbs, default_slurm
  from .. import jobs
  from . import _get_current_job_params, saveto
  ip = self.api
  ip.user_ns.pop("_lada_error", None)

  # creates mppalloc function.
  try: mppalloc = ip.ev(event.nbprocs)
  except Exception as e: 
    print "Could not make sense of --nbprocs argument {0}.\n{1}".format(event.nbprocs, e)
    return
  if mppalloc == None:
    def mppalloc(job): 
      """ Returns number of processes for this job. """
      N = len(job.structure.atoms) # number of atoms.
      if N % 2 == 1: N -= 1
      return N  

  # gets walltime.
  if re.match("\s*(\d{1,3}):(\d{1,2}):(\d{1,2})\s*", event.walltime) == None:
    try: walltime = ip.ev(event.walltime)
    except Exception as e: 
      print "Could not make sense of --walltime argument {0}.\n{1}".format(event.walltime, e)
      return
  else: walltime = event.walltime
  walltime = re.match("\s*(\d{1,3}):(\d{1,2}):(\d{1,2})\s*", walltime)
  if walltime != None:
    a, b, c = walltime.group(1), walltime.group(2), walltime.group(3)
    walltime = "{0:0>2}:{1:0>2}:{2:0>2}".format(a, b, c)
  else: 
    print "Could not make sense of --walltime argument {0}.".format(event.walltime)
    return

  # gets queue (aka partition in slurm), if any.
  kwargs = {}
  if event.__dict__.get('queue', None) != None: kwargs['queue'] = event.queue

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

  # where are we? important for default template.
  which = "SNLCLUSTER" in environ
  if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
  template = default_slurm if which else default_pbs
  # gets python script to launch in pbs.
  pyscript = jobs.__file__.replace(splitpath(jobs.__file__)[1], "runone.py")

  pbsscripts = []
  for current, path in jobdicts:
    # creates directory.
    with Changedir(path + ".pbs") as pwd: pass 
    # creates pbs scripts.
    for i, (name, job) in enumerate(current.root.iteritems()):
      if job.is_tagged: continue
      # does not launch successful jobs.
      if hasattr(job.functional, 'Extract') and not event.force: 
        p = join(dirname(path), name)
        extract = job.functional.Extract(p)
        if extract.success:
          print "Job {0} completed successfully. It will not be relaunched.".format(name)
          continue
      mppwidth = mppalloc(job) if hasattr(mppalloc, "__call__") else mppalloc
      if len(name) == 0: name = "{0}-root".format(splitpath(path)[1])
      name = name.replace("/", ".")
      pbsscripts.append(join(path+".pbs", name + ".pbs"))
      with open(pbsscripts[-1], "w") as file: 
        template( file, outdir=abspath(splitpath(path)[0]), jobid=i, mppwidth=mppwidth, name=name,\
                  pickle=splitpath(path)[1], pyscript=pyscript, ppath=".", walltime=walltime,\
                  **kwargs )
    print "Created scattered jobs in {0}.pbs.".format(path)

  if event.nolaunch: return
  # otherwise, launch.
  for script in pbsscripts:
    ip.system("{0} {1}".format("sbatch" if which else "qsub", script))



def launch(self, event):
  """ Launches PBS/slurm jobs.


      The usual setup is to first explore a dictionary of some sort, 
      then modify it and save it, and finally launch it.

      >>> # to recompute errors.
      >>> %explore errors path/to/original/pickle
      >>> %goto next
      >>> # modify dictionary here.
      >>> %showme functional
      >>> ...
      >>> %goto next
      >>> %showme functional
      >>> ...
      >>> # Saves the modified job.
      >>> # the new path could a different filename in the same directory.
      >>> # This way, the unsuccessful output from the first run will be
      >>> # overwritten.
      >>> %savejobs path/to/original/pickle.errors
      >>> # then  launch.
      >>> %launch scattered --walltime "24:00:00"
  """ 

  import argparse
  from os import environ
  from .. import queues as lada_queues

  which = "SNLCLUSTER" in environ
  if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
  queue = "--account" if which else "--queue"

  # main parser
  parser = argparse.ArgumentParser(prog='%launch')
  # options supported by all.
  opalls = argparse.ArgumentParser(add_help=False)
  opalls.add_argument('--walltime', type=str, default="05:59:59", \
                      help='walltime for jobs. Should be in hh:mm:ss format.')
  opalls.add_argument( 'pickle', metavar='FILE', type=str, nargs='*', default="", 
                       help='Optional path to a jobdictionary. If not present, the '\
                            'currently loaded job-dictionary will be launched.')


  # subparsers
  subparsers = parser.add_subparsers(help='Launches one job per untagged calculations')

  # launch scattered.
  scattered = subparsers.add_parser('scattered', description='Each calculation is a separate job.',\
                                    parents=[opalls])
  scattered.add_argument( '--nbprocs', type=str, default="None", dest="nbprocs",
                          help="Can be an integer, in which case it specifies "\
                               " the number of processes to exectute jobs with. "\
                               "Can also be a callable taking a JobDict as " \
                               "argument and returning a integer. Will default "\
                               "to as many procs as there are atoms in that particular structure.")
  scattered.add_argument('--nolaunch', action="store_true", dest="nolaunch")
  scattered.add_argument('--force', action="store_true", dest="force", \
                         help="Launches all untagged jobs, even those which completed successfully.")
  if which:
    scattered.add_argument( '--account', dest="account", choices=lada_queues,
                            default=(lada_queues[0] if len(lada_queues) > 0 else 'BES000'),
                            help="Account on which to launch job. Defaults to system default." )
  elif len(lada_queues) != 0: 
    scattered.add_argument( '--queue', dest="queue", choices=lada_queues, default=lada_queues[0],
                            help="Queue on which to launch job. Defaults to system default." )
  else: scattered.add_argument('--queue', dest="queue", type=str)
  scattered.set_defaults(func=launch_scattered)


  # parse arguments
  try: args = parser.parse_args(event.split())
  except SystemExit as e: return None
  else: args.func(self,args)


def launch_completer(self, info):
  """ Completion for launchers. """
  from os import environ
  from ._explore import _glob_job_pickles
  from .. import queues as lada_queues

  which = "SNLCLUSTER" in environ
  if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
  queue = "--account" if which else "--queue"
  
  ip = self.api
  data = info.line.split()
  if len(data)  <= 2 and data[-1] not in ["scattered"]: return ["scattered"] 
  if data[1] == "scattered": 
    if    (len(info.symbol) == 0 and data[-1] == "--walltime") \
       or (len(info.symbol) > 0  and data[-2] == "--walltime"):
      return [u for u in ip.user_ns if u[0] != '_' and isinstance(ip.user_ns[u], str)]
    if    (len(info.symbol) == 0 and data[-1] == "--nbprocs") \
       or (len(info.symbol) > 0  and data[-2] == "--nbprocs"):
      result = [u for u in ip.user_ns if u[0] != '_' and isinstance(ip.user_ns[u], int)]
      result.extend([u for u in ip.user_ns if u[0] != '_' and hasattr(u, "__call__")])
      return result
    if    (len(info.symbol) == 0 and data[-1] == queue) \
       or (len(info.symbol) > 0  and data[-2] == queue):
      return lada_queues
    result = ['--force', '--walltime', '--nbprocs', '--help']
    if len(lada_queues) > 0: result.append(queue) 
    result = list(set(result) - set(data))
    result.extend( _glob_job_pickles(ip, info.symbol) )
    return result
  raise IPython.ipapi.TryNext
         
