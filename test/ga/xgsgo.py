import random
from lada.ga.stepped import StepGA
from lada.ga.xgsgo import cut_and_splice, mix_atoms, mix_poscars, \
                          rand_struct_gen, jiggle_structure, \
                          update_fitness as update_fitness_xgsgo, \
                          objective as objective_xgsgo

def Extract(outdir=None, comm=None):
  from os.path import exists
  from os import getcwd
  from collections import namedtuple
  from pickle import load
  from numpy import array
  from lada.opt import Changedir

  if outdir == None: outdir = getcwd()
  Extract = namedtuple('Extract', ['success', 'directory', 'structure',
                                   'stoichiometry', 'total_energy'])
  if not exists(outdir): return Extract(False, outdir, None, None, None)
  with Changedir(outdir) as pwd:
    if not exists('OUTCAR'): return Extract(False, outdir, None, None, None)
    with open('OUTCAR', 'r') as file: indiv, value = load(file)
  types = set([a.type for a in indiv.atoms])
  stoic = array([len([0 for a in indiv.atoms if a.type == t]) for t in types])
  return Extract(True, outdir, indiv, stoic, value)

def functional(structure, outdir, comm, external=None):
  from pickle import dump
  from quantities import eV
  from lada.opt import Changedir

  with Changedir(outdir) as pwd:
    with open('OUTCAR', 'w') as file: dump((structure, (random.random()*100.0-60.0)*eV), file)

  return Extract(outdir)
functional.Extract = Extract


class XGSGO(StepGA):
  """ X-GSGO genetic algorithm metafunctional. """
  def __init__(self, functional, species, directory=None, 
               cns_rate=1e0, mix_atoms_rate=1e0, mix_poscar_rate=1e0,
               mutation_rate=0.5 ): 
    super(XGSGO, self).__init__(directory)
    self.matingops.add( cut_and_splice, rate=cns_rate )
    self.matingops.add( mix_atoms, rate=mix_atoms_rate )
    self.matingops.add( mix_poscars, rate=mix_poscar_rate )
    self.matingops.add( jiggle_structure, rate=mutation_rate )

    self.functional = functional
    """ Functional with which to get total energy. """
    self.species = list(species)
    """ Species in GA optimization. """
    # save it again.
    self._update_func()
    
  def random_individual(self, **kwargs):
    """ Generates a random individual. """
    return rand_struct_gen(species=self.species, **kwargs)

  update_fitness = update_fitness_xgsgo
  """ Finds convex-hull and figures out actual fitnesses.
  
      Updates fitness of all individuals.
  """

  def jobinator(self, job, indiv):
    """ Creates a job to compute an individual. """
    job.functional = self.functional
    job.jobparams['structure'] = indiv

  objective = objective_xgsgo
  """ Stores in individual information needed for convex-hull. 

      Does not actuall compute fitness, since we need all individuals for that.
  """


def random_indivs(n=20, N=20):
  from random import randint, random
  from IPython.ipapi import get
  from lada.ga.stepped import load

  ip = get()
  if 'ga' not in ip.user_ns: load()
  ga = ip.user_ns['ga']

  for i in xrange(n):
    M = randint(2, N)
    na = int(random()*M)
    print na, M-na
    ga.new_random_individuals(1, stoichiometry=[na, M-na])

