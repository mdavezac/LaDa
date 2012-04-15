""" Launches scattered calculations.

 
    This launch strategy will send one pbs/slurm job per lada job.

    >>> %launch scattered --walltime 24:00:00 
"""
__docformat__ = "restructuredtext en"

def launch(self, event, jobfolders):
  """ Launch scattered jobs: one job = one pbs script. """
  import re
  from copy import deepcopy
  from os.path import split as splitpath, join, dirname
  from ...misc import Changedir
  from ...jobs import __file__ as jobs_filename
  from ... import pbs_string, default_pbs, debug_queue, qsub_exe, default_comm

  comm = deepcopy(default_comm)
  comm['n'] = event.nbprocs
  comm["ppn"] = event.ppn
  
  pbsargs = deepcopy(default_pbs)
  pbsargs['comm'] = comm
  # creates mppalloc function.
  try: mppalloc = self.ev(event.nbprocs)
  except Exception as e: 
    print "Could not make sense of --nbprocs argument {0}.\n{1}".format(event.nbprocs, e)
    return
  if mppalloc is None:
    def mppalloc(job): 
      """ Returns number of processes for this job. """
      N = len(job.structure) # number of atoms.
      if N % 2 == 1: N -= 1
      return max(N, 1)  

  # gets walltime.
  if re.match("\s*(\d{1,3}):(\d{1,2}):(\d{1,2})\s*", event.walltime) is None:
    try: walltime = ip.ev(event.walltime)
    except Exception as e: 
      print "Could not make sense of --walltime argument {0}.\n{1}".format(event.walltime, e)
      return
  else: walltime = event.walltime
  walltime = re.match("\s*(\d{1,3}):(\d{1,2}):(\d{1,2})\s*", walltime)
  if walltime is not None:
    a, b, c = walltime.group(1), walltime.group(2), walltime.group(3)
    walltime = "{0:0>2}:{1:0>2}:{2:0>2}".format(a, b, c)
  else: 
    print "Could not make sense of --walltime argument {0}.".format(event.walltime)
    return

  # gets queue (aka partition in slurm), if any.
  pbsargs.update(comm)
  if event.__dict__.get("queue", None) is not None: pbsargs["queue"] = event.queue
  if event.__dict__.get("account", None) is not None: pbsargs["account"] = event.account
  if event.debug:
    assert debug_queue is not None, RuntimeError("debug_queue global variable has not been set.")
    pbsargs[debug_queue[0]] = debug_queue[1]

  # gets python script to launch in pbs.
  pyscript = jobs_filename.replace(splitpath(jobs_filename)[1], "runone.py")

  pbsscripts = []
  for current, path in jobfolders:
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
      pbsargs['n'] = mppalloc(job) if hasattr(mppalloc, "__call__") else mppalloc
      pbsargs['nnodes'] = (pbsargs['n'] + pbsargs['ppn'] - 1) // pbsargs['n']
      if len(name) == 0: pbsargs['name'] = "{0}-root".format(splitpath(path)[1])
      pbsargs['name'] = name.replace("/", ".")
      pbsscripts.append(join(path+".pbs", pbsargs['name'] + ".pbs"))
      if getattr(event, "prefix", None): pbsargs['name'] = "{0}-{1}".format(event.prefix, name)
      pbsargs['err'] = join("{0}.pbs".format(path), "err.{0}.pbs".format(pbsargs['name']))
      pbsargs['out'] = join("{0}.pbs".format(path), "out.{0}.pbs".format(pbsargs['name']))
      pbsargs['directory'] = dirname(path)
      pbsargs['scriptcommand'] = "{0} --nbprocs {n} --ppn {ppn} --jobid={1} {2}"\
                                 .format(pyscript, name, path, **pbsargs)
      with open(pbsscripts[-1], "w") as file: file.write(pbs_string.format(**pbsargs))
    print "Created {0} scattered jobs in {1}.pbs.".format(i, path)

  if event.nolaunch: return
  # otherwise, launch.
  for script in pbsscripts:
    self.system("{0} {1}".format(qsub_exe, script))


def completer(self, info, data):
  """ Completer for scattered launcher. """
  from ... import queues, accounts, debug_queue
  from .. import jobfolder_file_completer
  from ... import accounts, queues
  if len(data) > 0: 
    if data[-1] == "--walltime":
      return [u for u in self.user_ns if u[0] != '_' and isinstance(self.user_ns[u], str)]
    elif data[-1] == "--nbprocs": 
      result = [u for u in self.user_ns if u[0] != '_' and isinstance(self.user_ns[u], int)]
      result.extend([u for u in self.user_ns if u[0] != '_' and hasattr(u, "__call__")])
      return result
    elif data[-1] == '--ppn': return ['']
    elif data[-1] == "--prefix": return ['']
    elif data[-1] == "--queue": return queues
    elif data[-1] == "--account": return accounts
  result = ['--force', '--walltime', '--nbprocs', '--help']
  if len(queues) > 0: result.append("--queue") 
  if len(accounts) > 0: result.append("--account") 
  if debug_queue is not None: result.append("--debug")
  result.extend(jobfolder_file_completer(self, [info.symbol]))
  result = list(set(result) - set(data))
  return result

def parser(self, subparsers, opalls):
  """ Adds subparser for scattered. """ 
  from ... import queues, accounts, debug_queue, default_pbs, default_comm, qsub_exe
  result = subparsers.add_parser( 'scattered', 
                                  description="A separate PBS/slurm script is created for each "\
                                              "and every calculation in the job-folder "\
                                              "(or dictionaries).",
                                  parents=[opalls])
  result.add_argument('--walltime', type=str, default=default_pbs['walltime'], \
                         help='walltime for jobs. Should be in hh:mm:ss format. '\
                              'Defaults to ' + default_pbs['walltime'] + '.')
  result.add_argument('--prefix', action="store", type=str, help="Adds prefix to job name.")
  result.add_argument( '--nolaunch', action="store_true", dest="nolaunch",
                       help='Does everything except calling {0}.'.format(qsub_exe) )
  result.add_argument( '--nbprocs', type=str, default="None", dest="nbprocs",
                       help="Can be an integer, in which case it specifies "\
                            "the number of processes to exectute jobs with. "\
                            "Can also be a callable taking a JobFolder as " \
                            "argument and returning a integer. Will default "\
                            "to as many procs as there are atoms in that particular structure. "\
                            "Defaults to something close to the number of atoms in "
                            "the structure (eg good for VASP). ")
  result.add_argument( '--ppn', dest="ppn", default=default_comm.get('ppn', 1), type=int,
                       help="Number of processes per node. Defaults to {0}."\
                            .format(default_comm.get('ppn', 1)))
  if len(accounts) != 0:
    result.add_argument( '--account', dest="account", choices=accounts, default=accounts[0],
                         help="Account on which to launch job. Defaults to {0}.".format(accounts[0]) )
  else: result.add_argument('--account', dest="account", type=str,
                            help="Launches jobs on specific account if present.")
  if len(queues) != 0: 
    result.add_argument( '--queue', dest="queue", choices=queues, default=queues[0],
                         help="Queue on which to launch job. Defaults to {0}.".format(queues[0]) )
  else: result.add_argument('--queue', dest="queue", type=str,
                            help="Launches jobs on specific queue if present.")
  if debug_queue is not None:
    result.add_argument( '--debug', dest="debug", action="store_true", 
                         help="launches in interactive queue if present." )
  result.set_defaults(func=launch)
  return result
