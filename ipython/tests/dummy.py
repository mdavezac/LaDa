import random
import uuid

random.seed(12212)
def Extract(outdir=None):
  from os.path import exists
  from os import getcwd
  from collections import namedtuple
  from pickle import load
  from pylada.misc import Changedir

  if outdir == None: outdir = getcwd()
  Extract = namedtuple('Extract', ['success', 'directory', 'indiv', 'functional'])
  if not exists(outdir): return Extract(False, outdir, None, functional)
  with Changedir(outdir) as pwd:
    if not exists('OUTCAR'): return Extract(False, outdir, None, functional)
    with open('OUTCAR', 'r') as file: indiv, value = load(file)
  return Extract(True, outdir, indiv, functional)
  

def functional(indiv, outdir=None, value=False, **kwargs):
  from pylada.misc import Changedir
  from pickle import dump

  with Changedir(outdir) as pwd:
    with open('OUTCAR', 'w') as file: dump((indiv, value), file)

  return Extract(outdir)
functional.Extract = Extract


