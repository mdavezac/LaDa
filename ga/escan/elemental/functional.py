""" GA as a functional, for use with ladajobs. """
__docformat__ = "restructuredtext en"
from .extract import Extract as GAExtract

__all__ = ['Darwin']

class Darwin: 
  Extract = GAExtract
  """ Holds all GA parameters """
  RESTARTCAR = "GA_RESTARTCAR" 
  """ Pickle for restarting GA. """
  SAVECAR = "GA_SAVECAR" 
  """ Pickle to which GA saves. """
  def __init__(self, evaluator, **kwargs): 
    """ Initializes a GA functional. 
         
        :Parameters:
          evaluator : `lada.ga.escan.elemental.Bandgap`
            Functional which uses vff and escan for evaluating physical properties.
        :Kwarg crossover:
          Crossover object. Defaults to `lada.ga.escan.elemental.Crossover` ``(rate=0.5)``.
        :Kwarg mutation:
          Mutation object. `lada.ga.elemental.Mutation` ``(rate=2e0/float(len(.evaluator))``
        :Kwarg cm_rate:
          Crossover vs Mutation rate. Defaults to 0.1.
        :Kwarg rate: 
          Offspring rate. Defaults to ``max(0.1, 1/popsize)``.
        :Kwarg popsize: 
          Population size. Defaults to 100.
        :Kwarg maxgen: 
          Maximum number of generations. Defaults to 0.
        :Kwarg mean_conc:
          Mean concentration of individuals on random initialization. 
        :Kwarg stddev_conc:
          Standard deviation of concentration of individuals on random initialization.
        :Kwarg current_gen:
          Current generation. Defaults to 0 or to whatever is in the restart.
    """
    from ...standard import print_offspring, average_fitness, best
    from .. import CompareSuperCells
    from . import Crossover, Mutation
    from . import Individual

    self.checkpoints = [ print_offspring, 
                         average_fitness,
                         best,
                         self.__class__.print_nb_evals,
                         self.__class__.save]
    """ Checkpoints functions. """

    self.evaluator = evaluator
    """ Evaluator object taking an individual and computing its fitness. """
    Individual.size = len(self.evaluator)
    lattice = self.evaluator.escan.vff.lattice
    unit_cell = self.evaluator.converter.structure.cell
    self.compare = CompareSuperCells(lattice, unit_cell, self.evaluator.converter)
    """ Compares two different individuals. """


    # mating operations
    self.cm_rate = kwargs.pop("cm_rate", 0.1)
    """ Crossover versus Mutation rate. """
    self.crossover = kwargs.pop("crossover", Crossover(rate=0.5))
    """ Crossover operator. """
    self.mutation = kwargs.pop("mutation", Mutation(rate=2e0/float(len(self.evaluator))))
    """ Mutation operator. """

    # other parameters
    self.popsize = kwargs.pop("popsize", 100)
    self.rate    = kwargs.pop("rate", max(1e0/(float(self.popsize)-0.5), 0.1))
    self.max_gen = kwargs.pop("maxgen", 0)

    # Parameters for starting population.
    self.mean_conc = kwargs.pop("mean_conc", 0.5)
    """ Mean concentration of individuals on random initialization. """
    self.stddev_conc = kwargs.pop("stddev_conc", 0.5)
    """ Standard deviation of concentration of individuals on random initialization. """

    self.current_gen = kwargs.pop("current_gen", 0)
    """ Current generation. """

    self.population = []
    """ Current population """
    self.offspring = []
    """ Current offspring """

  def evaluation(self):
    """ Evaluates population. """
    from ...standard import population_evaluation
    population_evaluation(self, self.evaluator, self.comm, self.pools)

  def mating(self, *args, **kwargs):
    """ Calls mating operations. Removes hash code. """
    from .. import CompareSuperCells
    from ...standard import Mating
     
    # construct sequential mating operator.
    mating = Mating(sequential=False)
    mating.add(self.crossover, rate=cm_rate)
    mating.add(self.mutation,  rate=1-cm_rate) 
    # calls mating operator.
    result = mating(*args, **kwargs)
    # removes any has-code.
    CompareSuperCells.remove_hashcode(result)
    return result


  def taboo(self, indiv):
    """ No two individuals in the population and the offspring are the same. """
    from itertools import chain
    for a in chain(self.population, self.offspring):
      if self.compare(a, indiv): return True
    return False

  def Individual(self):
    """ Generates an individual (with mpi) """
    from numpy.random import normal, shuffle
    result = escan.elemental.Individual()
    if self.mean_conc == None or self.stddev_conc == None: return result
    N = 0
    while N == result.genes.size or N == 0:
      x = normal(self.mean_conc, self.stddev_conc)
      if x >= 1 or x <= 0e0: continue
      N = int(x * float(result.genes.size))
    result.genes[:N] = 1 # second type of appropriate site in lattice.
    result.genes[N:] = 0 # first type of appropriate site in lattice.
    shuffle(result.genes)
    return result

  @staticmethod
  def cmp_indiv(a, b, tolerance = 1e-12 ):
    """ Compares two individuals. 

        Maximizes fitness by default.
    """
    from math import fabs
    if fabs(a.fitness - b.fitness) <  tolerance: return 0
    return 1 if a.fitness < b.fitness else -1

  def restart(self):
    """ Saves current status. """
    import cPickle
    from os.path import exists
    from boost.mpi import broadcast
    if not exists(self.GA_RESTARTCAR):
      if self.comm.do_print == True:
        print "Restart file %s not found. Will start from scratch." % (self.GA_RESTARTCAR)
      return
    if self.comm.rank == 0:
      with open(self.GA_RESTARTCAR, "r") as file: this = cPickle.load(file)
      broadcast(self.comm, this,0)
      self.__dict__.update(this.__dict__)
    else:
      this = broadcast(self.comm, None,0)
      self.__dict__.update(this.__dict__)
    if self.comm.do_print:
      print "Restarting from file %s." % (self.GA_RESTARTCAR)
      print "Restarting at generation %i." % (self.current_gen)
      if len(self.population) == 0: print "Restarting with new population.\n"
      else:
        print "Restarting with population:"
        for indiv in self.population: print indiv, indiv.fitness

      if len(self.offspring) == 0: print "Restarting with empty offspring.\n"
      else:
        print "Restarting with offspring:"
        for indiv in self.offspring: print indiv, indiv.fitness

  def save(self):
    """ Saves current status. """
    import cPickle
    with open(self.GA_SAVECAR, "w") as file: cPickle.dump(self, file)
    return True
    
  def print_nb_evals(self):
    if self.comm.do_print: 
      print "Number of functional evaluations: ", self.evaluator.nbcalc

  def __call__(self, comm = None, **kwargs):
    from copy import deepcopy
    from boost.mpi import world
    from lada.ga import darwin as search
    # takes care of keyword arguments:
    if len(kwargs.keys()) > 0: 
      this = deepcopy(self)
      for key, value in kwargs:
        assert hasattr(self, key), TypeError("Unknown argument {0}.".format(key))
        setattr(self, key, value)
      this(comm=comm)
      return

    # mpi stuff
    self.comm = comm if comm != None else world
    self.comm.do_print = self.comm.rank == 0

    this.restart()
    search.run(this)
    del self.comm

  def __str__(self):
    return "Replacement rate: %f\nPopulation size: %i\nMaximum number of generations: %i\n" \
           "Size of each individuals: %i\n"\
           % (self.rate, self.popsize, self.max_gen, len(self.evaluator))

  def __getstate__(self):
    """ Returns current state. """
    d = self.__dict__.copy()
    d.pop("comm", None)
    return d

  def __setstate__(self, state):
    """ Resets current state. """
    self.__dict__.update(state)

  def __repr__(self):
    """ Returns representation of this instance. """
    return "from {0} import {1}\n{2}\n{3}\n{4}\n\n"\
           "functional = {1}(evaluator, crossover=crossover, mutation=mutation)\n"\
           "functional.cm_rate     = {5.cm_rate}\n"\
           "functional.popsize     = {5.popsize}\n"\
           "functional.max_gen     = {5.max_gen}\n"\
           "functional.current_gen = {5.current_gen}\n"\
           "functional.mean_conc   = {5.mean_conc}\n"\
           "functional.stddev_conc = {5.stddev_conc}\n"\
           .format( self.__class__.__module__, 
                    self.__class__.__name__,
                    repr(self.evaluator),
                    repr(self.crossover).replace("ga_operator", "crossover"),
                    repr(self.mutation).replace("ga_operator", "mutation"),
                    self )



        



