from ..functional import Functional as GAFunctional
class Functional(GAFunctional):
  """ XGSGO GA functional. """
  def __init__(self, functional, species, natoms=None, anions=None,
               cns_rate=1e0, mix_atoms_rate=1e0, mix_poscar_rate=1e0,
               mutation_rate=0.5, **kwargs): 
    super(XGSGO, self).__init__(directory)
    self.matingops.add( cut_and_splice, rate=cns_rate )
    self.matingops.add( mix_atoms, rate=mix_atoms_rate )
    self.matingops.add( mix_poscars, rate=mix_poscar_rate )
    self.matingops.add( jiggle_structure, rate=mutation_rate )

    self.functional = functional
    """ Functional with which to get total energy. """
    self.species    = list(species)
    """ Species in GA optimization. """
    self.anions     = anions
    """ Anionic species in the calculation. """
    self.natom      = (2, 20) if natoms is None else (min(natoms), max(natoms))
    """ Maximum number of atoms. """
    if self.natoms[0] < 1: self.natoms = (1, self.natoms[1])

  def random_individual(self):
    """ Returns a new random individual. """
    from random import randint, shuffle
    from .random import random_structure, populate, populate_anion_cation
    from numpy import dot, array
    N = randint(*self.natoms)
    result = random_structure(N)
    natoms = []
    for i in xrange(len(self.species)):
      natoms.append(0 if natoms == 0 else randint(1, N))
      N -= natoms[-1]
    suffle(natoms)
    species = {}
    for s, n in zip(self.species, natoms): species[s] = N
    if self.anions is not None and len(self.anions) > 0:
      populate_anion_cation(structure, species, self.anions)
    else: populate(structure, species)
    return structure

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
    from .random import taboo
    if len(structure) < self.natoms[0] or len(structure) > self.natoms[1]: 
      return False
    return taboo(individual)

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

