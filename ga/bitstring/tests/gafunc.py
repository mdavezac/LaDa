from lada.ga.functional import Functional as GAFunctional
class Functional(GAFunctional):
  """ GA bitstring functional. """
  def __init__(self, program, size=10, maxgen=-1, **kwargs):
    """ Initializes functional. """
    from lada.ga import bitstring
    super(Functional, self).__init__(**kwargs)
    self.size = size
    self.program = program
    self.matingops.add( bitstring.Crossover(rate=0.25), rate=0.8 )
    self.matingops.add( bitstring.Mutation(rate=3e0/float(bitstring.Individual.size)), rate=0.2 )
    self.maxgen = maxgen

  def random_individual(self):
    """ Returns a new random individual. """
    from lada.ga.bitstring import Individual
    return Individual(size=self.size)

  def jobinator(self, job, indiv):
    """ Initializes a new job from an individual. """
    from functional import SerialFunctional
    from random import random
    job.functional = SerialFunctional(self.program, order=sum(indiv.genes))
    if random() > 0.8: # generate artificial failure
      job.functional.order = 666
  def objective(self, extract, indiv):
    """ Returns error. """
    return extract.error

  def checkpoints(self):
    """ Returns true when result is in population. """
    from numpy import all
    if self.generation > self.maxgen and self.maxgen >= 0: return True
    for indiv in self.population:
      if all(indiv.genes == 1): return True
    return False
  def is_taboo(self, individual, mechanics='mating'): return False

  @property 
  def process(self):
    """ Current process. 

        The process is an object which manages launching actual calculation.
        It basically takes a job-folder and launches calculations as
        appropriate. 

        By default, it creates a
        :py:class:`~lada.process.jobfolder.JobFolderProcess`. It might be a
        good idea to use something more specific, such as a
        :py:class:`~lada.process.pool.PoolProcess` instantiated with the
        correct ressource allocation function.

        The process is not saved when self is pickled. It is created anew each
        time this functional runs.
    """ 
    from lada import default_comm
    from lada.process import JobFolderProcess
    if '_process' not in self.__dict__: 
      self._process = JobFolderProcess(self.jobfolder, self.calcdir, nbpools=4)
    return self._process
