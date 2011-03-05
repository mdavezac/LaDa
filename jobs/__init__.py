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
           'default_params', 'JobParams']

from ..opt.decorators import broadcast_result
from .bleeder import Bleeder
from .jobdict import JobDict
from .manipulator import JobParams
from .extract import AbstractMassExtract, MassExtract, AbstractMassExtractDirectories
import default_params
DefaultParams = default_params

@broadcast_result(key=True)
def save(jobdict, path = None, overwrite=False): 
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
      :type comm: boost.mpi.comm
      :keyword overwrite: if True, then overwrites file.

      This method first acquire an exclusive lock on the file before writing
      (see `lada.opt.open_exclusive`).  This way not two processes can
      read/write to this file while using this function.
  """ 
  from os.path import exists
  from pickle import dump
  from ..opt import open_exclusive, RelativeDirectory
  if path == None: path = "pickled_jobdict"
  path = "pickled_jobdict" if path == None else RelativeDirectory(path).path
  if exists(path) and not overwrite: 
    print path, "exists. Please delete first if you want to save the job dictionary."
    return
  with open_exclusive(path, "wb") as file: dump(jobdict, file)
  print "Saved job dictionary to %s." % (path)

@broadcast_result(key=True)
def load(path = None): 
  """ Unpickles a job from file. 
 
      :keyword path: Filename of a pickled jobdictionary.
      :keyword comm: MPI processes for which to read job-dictionary.
      :type comm: boost.mpi.communicator
      :return: Returns a JobDict object.

      This method first acquire an exclusive lock (using os dependent lockf) on
      the file before reading. This way not two processes can read/write to
      this file while using this function.
  """ 
  from os.path import exists
  from pickle import load as load_pickle
  from ..opt import open_exclusive, RelativeDirectory
  path = "pickled_jobdict" if path == None else RelativeDirectory(path).path
  assert exists(path), IOError("File " + path + " does not exist.")
  with open_exclusive(path, "rb") as file: result = load_pickle(file)
  print "Loaded job list from", path, "."
  return result
