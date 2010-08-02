""" IPython launch magic function. """

def launch_scattered(self, event):
  """ Launch scattered jobs: one job = one pbs script. 
  
      A dictionary must have been loaded with explore, or saved. 
      Output will be saved in the current directory of the job.

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
      >>> %launch_scattered 

  """
  import re
  from os import environ
  from os.path import split as splitpath, join, exists
  from ..opt.changedir import Changedir
  from ..jobs.templates import default_pbs, default_slurm
  from .. import jobs
  from . import _get_current_job_params, saveto
  ip = self.api
  current, path = _get_current_job_params(self, 2)
  if current == None: return
  if path == None: return
  ip.user_ns.pop("_lada_error", None)

  # creates mppalloc function.
  if len(event.split()) == 0: 
    def mppalloc(job): 
      """ Returns number of processes for this job. """
      N = len(job.structure.atoms) # number of atoms.
      if N % 2 == 1: N -= 1
      return N  
  elif len(event.split()) > 1:
    print "launch_scattered accepts no more than one parameter."
    return
  elif event.split()[0] in ip.user_ns:
    mppalloc = ip.user_ns[event.split()[0]]
  else:
    try: mppalloc = int(event.split()[0])
    except: 
      print "Could not make sense of argument %s." % event
      return

  # where are we? important for default template.
  which = "SNLCLUSTER" in environ
  if which: which = environ["SNLCLUSTER"] in ["redrock", "redmesa"]
  template = default_slurm if which else default_pbs
  # gets python script to launch in pbs.
  pyscript = jobs.__file__.replace(splitpath(jobs.__file__)[1], "runone.py")
  # creates directory.
  with Changedir(path + ".pbs") as pwd: pass 
  # saving pickle
  saveto(self, path)
  if "_lada_error" in ip.user_ns:
    if ip.user_ns["_lada_error"] == "User said no save.":
      print "Job-dictionary not saved = jobs not launched."
    return
  # creates pbs scripts.
  pbsscripts = []
  for i, (job, name) in enumerate(current.walk_through()):
    if job.is_tagged: continue
    mppwidth = mppalloc(job) if hasattr(mppalloc, "__call__") else mppalloc
    name = name.replace("/", ".")
    pbsscripts.append(join(path+".pbs", name + ".pbs"))
    with open(pbsscripts[-1], "w") as file: 
      template( file, outdir=splitpath(path)[0], jobid=i, mppwidth=mppwidth, name=name,\
                pickle = path, pyscript=pyscript, ppath=splitpath(path)[0] )
  if event.find("nolaunch") != -1: return 
  # otherwise, launch.
  for script in pbsscripts:
    if which: ip.system("sbatch %s" % script)
    else: ip.system("qsub %s" % script)



