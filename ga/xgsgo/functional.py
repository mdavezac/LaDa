from ..functional import Functional as GAFunctional
from .objective import objective as objective_xgsgo
class Functional(GAFunctional):
  """ XGSGO GA functional. """
  def __init__(self, functional, species, natoms=None, anions=None,
               cns_rate=1e0, mix_atoms_rate=1e0, mix_poscar_rate=1e0,
               mutation_rate=0.5, **kwargs): 
    from .operators import cut_and_splice, mix_atoms, mix_poscars,\
                           jiggle_structure
    super(Functional, self).__init__()
    if cns_rate > 0: self.matingops.add(cut_and_splice, rate=cns_rate)
    if mix_atoms_rate > 0: self.matingops.add(mix_atoms, rate=mix_atoms_rate)
    if mix_poscar_rate > 0:
      self.matingops.add(mix_poscars, rate=mix_poscar_rate)
    if mutation_rate > 0:
      self.matingops.add(jiggle_structure, rate=mutation_rate)

    self.functional = functional
    """ Functional with which to get total energy. """
    self.species    = list(species)
    """ Species in GA optimization. """
    self.anions     = anions
    """ Anionic species in the calculation. """
    self.natoms     = (2, 20) if natoms is None else (min(natoms), max(natoms))
    """ Maximum number of atoms. """
    if self.natoms[0] < 1: self.natoms = (1, self.natoms[1])

  def random_individual(self):
    """ Returns a new random individual. """
    from random import randint, shuffle
    from numpy import dot, array
    from .initialization import random_structure, populate, populate_anion_cation
    N = randint(*self.natoms)
    result = random_structure(N)
    natoms = []
    for i in xrange(len(self.species)):
      natoms.append(N if N < 2 else randint(1, N))
      N -= natoms[-1]
    shuffle(natoms)
    species = {}
    for s, n in zip(self.species, natoms): species[s] = N
    if self.anions is not None and len(self.anions) > 0:
      populate_anion_cation(result, species, self.anions)
    else: populate(result, species)
    return result

  def jobinator(self, job, indiv):
    """ Initializes a new job from an individual. """
    job.functional = self.functional
    job.jobparams['structure'] = indiv
      
  objective = objective_xgsgo
  """ Stores in individual information needed for convex-hull. 

      Does not actuall compute fitness, since we need all individuals for that.
  """
  def is_taboo(self, individual):
    """ Checks whether an individual is taboo. """
    from .initialization import taboo
    if len(individual) < self.natoms[0] or len(individual) > self.natoms[1]: 
      return False
    return taboo(individual, same_first_neigh=-1)

  @property 
  def process(self):
    """ Current process. 

        The process is an object which manages launching actual calculation.
        It basically takes a job-folder and launches calculations as
        appropriate. 

        By default, it creates a
        :py:class:`~lada.process.pool.PoolProcess` where each job is allocated
        the N procs, where N is the even number closest from below to the
        number of atoms in the stucturel.

        The process is not saved when self is pickled. It is created anew each
        time this functional runs.
    """ 
    from lada.process import PoolProcess
    def nbprocs(job):
      return len(job.structure) - len(job.structure) % 1
    if '_process' not in self.__dict__: 
      self._process = PoolProcess(self.jobfolder, self.calcdir, nbprocs)
    return self._process

  def __call__(self, *args, **kwargs):
    """ Performs GA. """
    self.nbprocs = kwargs['comm']['n']
    return super(Functional, self).__call__(*args, **kwargs)

  def go_to_next_iteration(self):
    """ Returns True when going to next generation.

        Not all calculations need be finished.
        It is the job of this function to catch errors in process execution and
        correct them.
    """ 
    from time import sleep
    from ..process import Fail
    
    fewjobsleft = max(len(self.offspring) * 0.1, 0)
    lotsofcpusidle = min(self.process._alloc.values())
    while self.process.nbjobsleft > fewjobsleft \
          and self.process._comm['n'] < lotsofcpusidle:
      try:
        if self.process.poll(): break
      except Fail: break
      else: sleep(5)
    self._process_errors()
    return True
