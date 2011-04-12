""" Launches scattered calculations.

 
    This launch strategy will send one pbs/slurm job per lada job.

    >>> %launch scattered --walltime 24:00:00 
"""
__docformat__ = "restructuredtext en"

def launch(self, event, jobdicts):
  """ Launch scattered jobs: one job = one pbs script. """
  import re
  from os import environ
  from os.path import split as splitpath, join, exists, abspath, dirname
  from ...opt.changedir import Changedir
  from ...jobs.templates import default_pbs, default_slurm
  from ...jobs import __file__ as jobs_filename
  from .. import saveto
  ip = self.api
  ip.user_ns.pop("_lada_error", None)

  # where are we? important for default template.
  which = "SNLCLUSTER" in environ
  if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
  queue = "account" if which else "queue"
  template = default_slurm if which else default_pbs

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
  if event.__dict__.get(queue, None) != None: kwargs[queue] = getattr(event, queue)
  if which and event.debug: kwargs["partition"] = "inter"

  # gets python script to launch in pbs.
  pyscript = jobs_filename.replace(splitpath(jobs_filename)[1], "runone.py")

  if environ.get('NERSC_HOST', 'none') == 'hopper2':
    from os.path import basename, exists
    name = basename(environ['VIRTUAL_ENV'])
    pyscript = pyscript.replace('/{0}/'.format(name), '/{0}/'.format(name+"_mpi"))
    assert exists(pyscript),\
           RuntimeError("Could not find find file {0}.\n"\
                        "{1}_mpi environment does not exist?\n"\
                        "Or not actually on hopper2?".format(pyscript, name))
    kwargs['header'] = "\nsource {0}_mpi/bin/activate".format(environ['VIRTUAL_ENV'])

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
      if getattr(event, "prefix", None): name = "{0}-{1}".format(event.prefix, name)
      with open(pbsscripts[-1], "w") as file: 
        template( file, outdir=abspath(splitpath(path)[0]), jobid=i, mppwidth=mppwidth, name=name,\
                  pickle=splitpath(path)[1], pyscript=pyscript, ppath=".", walltime=walltime,\
                  **kwargs )
    print "Created scattered jobs in {0}.pbs.".format(path)

  if event.nolaunch: return
  # otherwise, launch.
  for script in pbsscripts:
    ip.system("{0} {1}".format("sbatch" if which else "qsub", script))


def completer(self, info, data, computer):
  """ Completer for scattered launcher. """
  from ... import queues as lada_queues
  from .._explore import _glob_job_pickles
  queue = "--account" if computer else "--queue"
  ip = self.api
  if    (len(info.symbol) == 0 and data[-1] == "--walltime") \
     or (len(info.symbol) > 0  and data[-2] == "--walltime"):
    return [u for u in ip.user_ns if u[0] != '_' and isinstance(ip.user_ns[u], str)]
  if    (len(info.symbol) == 0 and data[-1] == "--nbprocs") \
     or (len(info.symbol) > 0  and data[-2] == "--nbprocs"):
    result = [u for u in ip.user_ns if u[0] != '_' and isinstance(ip.user_ns[u], int)]
    result.extend([u for u in ip.user_ns if u[0] != '_' and hasattr(u, "__call__")])
    return result
  if    (len(info.symbol) == 0 and data[-1] == "--prefix") \
     or (len(info.symbol) > 0  and data[-2] == "--prefix"):
    return []
  if    (len(info.symbol) == 0 and data[-1] == queue) \
     or (len(info.symbol) > 0  and data[-2] == queue):
    return lada_queues
  result = ['--force', '--walltime', '--nbprocs', '--help']
  if len(lada_queues) > 0: result.append(queue) 
  if computer: result.append("--debug")
  result.extend( _glob_job_pickles(ip, info.symbol) )
  result = list(set(result) - set(data))
  return result

def parser(self, subparsers, which, opalls):
  """ Adds subparser for scattered. """ 
  from ... import queues as lada_queues
  result = subparsers.add_parser( 'scattered', 
                                  description="A separate PBS/slurm script is created for each "\
                                              "and every calculation in the jobdictionary "\
                                              "(or dictioanaries).",
                                  parents=[opalls])
  result.add_argument('--walltime', type=str, default="05:59:59", \
                         help='walltime for jobs. Should be in hh:mm:ss format.')
  result.add_argument( '--nbprocs', type=str, default="None", dest="nbprocs",
                       help="Can be an integer, in which case it specifies "\
                            " the number of processes to exectute jobs with. "\
                            "Can also be a callable taking a JobDict as " \
                            "argument and returning a integer. Will default "\
                            "to as many procs as there are atoms in that particular structure.")
  result.add_argument('--prefix', action="store", type=str, help="Adds prefix to job name.")
  result.add_argument('--nolaunch', action="store_true", dest="nolaunch")
  result.add_argument('--force', action="store_true", dest="force", \
                         help="Launches all untagged jobs, even those which completed successfully.")
  if which:
    result.add_argument( '--account', dest="account", choices=lada_queues,
                         default=(lada_queues[0] if len(lada_queues) > 0 else 'BES000'),
                         help="Account on which to launch job. Defaults to system default." )
    result.add_argument( '--debug', dest="debug", action="store_true", 
                         help="launches in interactive queue if present." )
  elif len(lada_queues) != 0: 
    result.add_argument( '--queue', dest="queue", choices=lada_queues, default=lada_queues[0],
                         help="Queue on which to launch job. Defaults to system default." )
  else: result.add_argument('--queue', dest="queue", type=str,
                            help="Launches jobs on specific queue if present.")
  result.set_defaults(func=launch)
  return result
