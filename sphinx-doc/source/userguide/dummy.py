def Extract(outdir=None):
  from os.path import exists
  from os import getcwd
  from collections import namedtuple
  from pickle import load
  from lada.misc import Changedir

  if outdir == None: outdir = getcwd()
  Extract = namedtuple('Extract', ['success', 'directory', 'structure', 'value', 'functional'])
  if not exists(outdir): return Extract(False, outdir, None, None, functional)
  with Changedir(outdir) as pwd:
    if not exists('OUTCAR'): return Extract(False, outdir, None, None, functional)
    with open('OUTCAR', 'r') as file: structure, value = load(file)
  return Extract(True, outdir, structure, value, functional)
  

def functional(structure, outdir=None, value=False, **kwargs):
  from lada.misc import Changedir
  from copy import deepcopy
  from pickle import dump

  structure = deepcopy(structure)
  structure.value = value
  with Changedir(outdir) as pwd:
    with open('OUTCAR', 'w') as file: dump((structure, value, functional), file)

  return Extract(outdir)
  return structure
functional.Extract = Extract


