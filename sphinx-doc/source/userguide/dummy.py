def Extract(outdir=None):
  """ An extraction function for a dummy functional """
  from os.path import exists
  from os import getcwd
  from collections import namedtuple
  from pickle import load
  from lada.misc import Changedir

  if outdir == None: outdir = getcwd()
  Extract = namedtuple('Extract', ['success', 'directory', 'energy', 'structure', 'value', 'functional'])
  if not exists(outdir): return Extract(False, outdir, None, None, None, None)
  with Changedir(outdir) as pwd:
    if not exists('OUTCAR'): return Extract(False, outdir, None, None, None, None)
    with open('OUTCAR', 'r') as file:
      structure, energy, value, functional = load(file)
      return Extract(True, outdir, energy, structure, value, functional)
  
def functional(structure, outdir=None, value=False, **kwargs):
  """ A dummy functional """
  from copy import deepcopy
  from pickle import dump
  from random import random
  from lada.misc import Changedir

  structure = deepcopy(structure)
  structure.value = value
  with Changedir(outdir) as pwd:
    with open('OUTCAR', 'w') as file: dump((random(), structure, value, functional), file)

  return Extract(outdir)
  return structure
functional.Extract = Extract

def create_jobs():
  """ Simple job-folders. """
  from lada.jobfolder import JobFolder
  from lada.crystal.binary import zinc_blende

  root = JobFolder()
  
  for name, value, species in zip( ['diamond', 'diamond/alloy', 'GaAs'], 
                                   [0, 1, 2],
                                   [('Si', 'Si'), ('Si', 'Ge'), ('Ga', 'As')] ):
    job = root / name 
    job.functional = functional
    job.params['value'] = value
    job.params['structure'] = zinc_blende()
    for atom, specie in zip(job.structure, species): atom.type = specie

  return root
