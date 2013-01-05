""" Tests launching jobs through magic function. """
from pylada.opt import AbstractExtractBase
from pylada.opt.decorators import broadcast_result

def functional(outdir=None, comm=None, **kwargs):
  """ Dummy functional. """
  from os import getcwd
  from os.path import join
  from pylada.opt import Changedir
  if outdir is None: outdir = getcwd()
  comm.barrier()
  with Changedir(outdir, comm) as cwd:
    if comm.is_root:
      with open('out', 'w') as file:
        file.write("comm size : {0}\n".format(comm.size))
        file.write("kargs : {0}\n".format(kwargs))

class Extract(AbstractExtractBase):
  """ Defines success for dummy functional. """
  @property
  @broadcast_result(attr=True, which=0)
  def success(self): 
    """ True if relevant file exists. """
    from os.path import exists, join
    return exists(join(self.directory, 'out'))

functional.Extract = Extract

def create_jobs(n=3):
  """ Returns dictionary with fake jobs. """
  from launch import functional
  from pylada.jobs import JobFolder

  root = JobFolder()
  for i in xrange(n):
    job = root / "results" / str(i)
    job.functional = functional
    job.jobparams["param"] = 'param'

  return root
