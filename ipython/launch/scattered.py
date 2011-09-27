""" Launches scattered calculations.

 
    This launch strategy will send one pbs/slurm job per lada job.

    >>> %launch scattered --walltime 24:00:00 
"""
__docformat__ = "restructuredtext en"

def launch(self, event, jobdicts):
  """ Launch scattered jobs: one job = one pbs script. """
  import re
  from os.path import split as splitpath, join, abspath, dirname
  from ...opt.changedir import Changedir
  from ...jobs import __file__ as jobs_filename
  from ... import template_pbs, debug_queue, qsub_exe
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
      return max(N, 1)  

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
  kwargs = { "ppernode": event.ppn,
             "memlim": event.memlim,
             "external": event.external } 
  if event.__dict__.get("queue", None) != None: kwargs["queue"] = getattr(event, "queue")
  if event.__dict__.get("account", None) != None: kwargs["account"] = getattr(event, "account")
  if event.debug:
    assert debug_queue != None, RuntimeError("debug_queue global variable has not been set.")
    kwargs[debug_queue[0]] = debug_queue[1]

  # gets python script to launch in pbs.
  pyscript = jobs_filename.replace(splitpath(jobs_filename)[1], "runone.py")

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
        template_pbs( file, outdir=abspath(splitpath(path)[0]), jobid=i, mppwidth=mppwidth, name=name,\
                      pickle=splitpath(path)[1], pyscript=pyscript, ppath=".", walltime=walltime,\
                      **kwargs )
    print "Created scattered jobs in {0}.pbs.".format(path)

  if event.nolaunch: return
  # otherwise, launch.
  for script in pbsscripts:
    ip.system("{0} {1}".format(qsub_exe, script))


def completer(self, info, data):
  """ Completer for scattered launcher. """
  from ... import queues, accounts, debug_queue
  from .._explore import _glob_job_pickles
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
  if    (len(info.symbol) == 0 and data[-1] == "--queue") \
     or (len(info.symbol) > 0  and data[-2] == "--queue"):
    return queues
  if    (len(info.symbol) == 0 and data[-1] == "--account") \
     or (len(info.symbol) > 0  and data[-2] == "--account"):
    return accounts
  result = ['--force', '--walltime', '--nbprocs', '--help', '--external']
  if len(queues) > 0: result.append("--queue") 
  if len(accounts) > 0: result.append("--account") 
  if debug_queue != None: result.append("--debug")
  result.extend( _glob_job_pickles(ip, info.symbol) )
  result = list(set(result) - set(data))
  return result

def parser(self, subparsers, opalls):
  """ Adds subparser for scattered. """ 
  from ... import queues, accounts, debug_queue, default_walltime
  from ... import cpus_per_node
  result = subparsers.add_parser( 'scattered', 
                                  description="A separate PBS/slurm script is created for each "\
                                              "and every calculation in the jobdictionary "\
                                              "(or dictioanaries).",
                                  parents=[opalls])
  result.add_argument('--walltime', type=str, default=default_walltime, \
                         help='walltime for jobs. Should be in hh:mm:ss format. '\
                              'Defaults to ' + default_walltime + '.')
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
  result.add_argument( '--memlim', dest="memlim", default="guess",
                       help="Memory limit per process imposed by ulimit. "\
                            "\"guess\" lets lada make an uneducated guess. ")
  result.add_argument( '--ppn', dest="ppn", default= cpus_per_node,
                       help="Number of processes per node.")
  if len(accounts) != 0:
    result.add_argument( '--account', dest="account", choices=accounts, default=accounts[0],
                         help="Account on which to launch job. Defaults to system default." )
  else: result.add_argument('--account', dest="account", type=str,
                            help="Launches jobs on specific account if present.")
  if len(queues) != 0: 
    result.add_argument( '--queue', dest="queue", choices=queues, default=queues[0],
                         help="Queue on which to launch job. Defaults to system default." )
  else: result.add_argument('--queue', dest="queue", type=str,
                            help="Launches jobs on specific queue if present.")
  if debug_queue != None:
    result.add_argument( '--debug', dest="debug", action="store_true", 
                         help="launches in interactive queue if present." )
  result.set_defaults(func=launch)
  return result
