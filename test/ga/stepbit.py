import random
import uuid
from lada.ga.stepped import StepGA


def uuid4(thislist=[]):
  if len(thislist) == 0:
    letters=list('abcdef0123456789')
    thislist.append(0)
    for i in xrange(1000):
      dummy = ''
      for j in xrange(32): dummy += random.choice(letters)
      thislist.append(uuid.UUID(dummy))
    
  thislist[0] += 1
  return thislist[thislist[0]]

random.seed(12212)
uuid.uuid4 = uuid4
uuid.uuid4()
def Extract(outdir=None, comm=None):
  from os.path import exists
  from os import getcwd
  from collections import namedtuple
  from pickle import load
  from lada.opt import Changedir

  if outdir == None: outdir = getcwd()
  Extract = namedtuple('Extract', ['success', 'directory', 'indiv'])
  if not exists(outdir): return Extract(False, outdir, None)
  with Changedir(outdir) as pwd:
    if not exists('OUTCAR'): return Extract(False, outdir, None)
    with open('OUTCAR', 'r') as file: indiv, value = load(file)
  return Extract(True, outdir, indiv)
  

def functional(indiv, outdir, comm=None, external=None):
  from lada.opt import Changedir
  from pickle import dump

  with Changedir(outdir) as pwd:
    with open('OUTCAR', 'w') as file: dump((indiv, random.random()), file)

  return Extract(outdir)
functional.Extract = Extract


class BitString(StepGA):
  def __init__(self, directory=None): 
    from lada.ga.bitstring import Crossover, Mutation
    super(BitString, self).__init__(directory)
    self.matingops.add( Crossover(rate=0.5), rate=1 )
    self.matingops.add( Mutation(rate=2e0/float(self.random_individual().size)), rate=0.2 )
    # save it again.
    self._update_func()
    
    
  def random_individual(self, **kwargs):
    from lada.ga.bitstring import Individual
    return Individual()
  def jobinator(self, job, indiv):
    job.functional = functional
    job.jobparams['indiv'] = indiv
  def objective(self, extract, indiv, _target = []):
    from numpy import array, sum, abs
    if len(_target) == 0:
      _target.append(array([random.choice([0,1]) for i in xrange(len(indiv.genes))]))
    return sum(abs(indiv.genes - _target[0]))
