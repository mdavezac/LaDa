""" Classes and functions pertaining to job-management. 

    Jobs are described within job-dictionaries. These dictionaries resemble
    directory trees, where some directories contain input files (eg job
    parameters) and are marked for execution.

    In addition, this package contains classes, such as MassExtract and
    JobParams, capable of navigating output and input from jobdictionaries.
"""
__docformat__ = "restructuredtext en"
__all__ = ['JobDict', 'walk_through', 'save', 'load', 'MassExtract',
           'AbstractMassExtract', 'AbstractMassExtractDirectories', 'Bleeder',
           'JobParams', 'SuperCall']

from .bleeder import Bleeder
from .jobdict import JobDict, SuperCall
from .manipulator import JobParams
from .extract import AbstractMassExtract, AbstractMassExtractDirectories
from .massextract import MassExtract 

def save(jobdict, path=None, overwrite=False, timeout=None): 
  """ Pickles a job to file. 
 
      :keyword jobdict: A job-dictionary to pickle. 
      :type jobdict: `JobDict`
      :keyword path: 
          filename of file to which to save pickle. overwritten. If None then
          saves to "pickled_jobdict"
      :type path: str or None
      :keyword comm:
        Convenience parameter. Only root process actually saves.
        Other processes wait silently.
      :type comm: `mpi.Communicator`
      :keyword overwrite: if True, then overwrites file.

      This method first acquire an exclusive lock on the file before writing
      (see `lada.opt.open_exclusive`).  This way not two processes can
      read/write to this file while using this function.
  """ 
  from os.path import exists
  from pickle import dump
  from ..misc import open_exclusive, RelativePath
  from .. import is_interactive
  if path is None: path = "pickled_jobdict"
  path = "job.dict" if path is None else RelativePath(path).path
  if exists(path) and not overwrite: 
    if is_interactive:
      print path, "exists. Please delete first if you want to save the job dictionary."
      return
    else: raise IOError('{0} already exists. By default, will not overwrite.'.format(path))
  with open_exclusive(path, "wb", timeout=None) as file: dump(jobdict, file)
  if is_interactive: print "Saved job dictionary to {0}.".format(path)

def load(path = None, timeout=None): 
  """ Unpickles a job from file. 
 
      :keyword path: Filename of a pickled jobdictionary.
      :keyword comm: MPI processes for which to read job-dictionary.
      :type comm: `mpi.Communicator`
      :return: Returns a JobDict object.

      This method first acquire an exclusive lock (using os dependent lockf) on
      the file before reading. This way not two processes can read/write to
      this file while using this function.
  """ 
  from os.path import exists
  from pickle import load as load_pickle
  from ..misc import open_exclusive, RelativePath
  from .. import is_interactive
  path = "job.dict" if path is None else RelativePath(path).path
  if not exists(path): raise IOError("File " + path + " does not exist.")
  with open_exclusive(path, "rb", timeout=timeout) as file: result = load_pickle(file)
  if is_interactive: print "Loaded job list from {0}.".format(path)
  return result
