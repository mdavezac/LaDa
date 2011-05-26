""" GA as a functional, for use with ladajobs. """
__docformat__ = "restructuredtext en"
from ..functional import Darwin as EscanDarwin

__all__ = ['Darwin']

class Darwin(EscanDarwin): 
  """ GA functional for optimizations of epitaxial structures. """

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
    from .. import CompareSuperCells
    from . import Crossover, Mutation
    super(Darwin, self).__init__(evaluator, **kwargs)

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

    # Parameters for starting population.
    self.mean_conc = kwargs.pop("mean_conc", 0.5)
    """ Mean concentration of individuals on random initialization. """
    self.stddev_conc = kwargs.pop("stddev_conc", 0.5)
    """ Standard deviation of concentration of individuals on random initialization. """

  def mating(self, *args, **kwargs):
    """ Calls mating operations. Removes hash code. """
    from .. import CompareSuperCells
    from ...standard import Mating
     
    # construct sequential mating operator.
    mating = Mating(sequential=False)
    mating.add(self.crossover, rate=self.cm_rate)
    mating.add(self.mutation,  rate=1-self.cm_rate) 
    # calls mating operator.
    result = mating(self, *args, **kwargs)
    # removes any hash-code.
    CompareSuperCells.remove_hashcode(result)
    return result


  def Individual(self):
    """ Generates an individual (with mpi) """
    from . import Individual
    from numpy.random import normal, shuffle
    result = Individual(size=len(self.evaluator))
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

  def __repr__(self):
    """ Returns representation of this instance. """
    header = "from {0.__class__.__module__} import {0.__class__.__name__}\n".format(self)
    string = repr(self.evaluator) + "\n" + repr(self.matingops) + "\n"
    string += "functional = {0.__class__.__name__}(evaluator)\n"\
              "functional.pools       = {0.pools}\n"\
              "functional.nmax        = {0.nmax}\n"\
              "functional.nmin        = {0.nmin}\n"\
              "functional.popsize     = {0.popsize}\n"\
              "functional.max_gen     = {0.max_gen}\n"\
              "functional.current_gen = {0.current_gen}\n"\
              "functional.rate        = {0.rate}\n"\
              "functional.rootworkdir = {1}\n"\
              "functional.current_gen = {0.current_gen}\n"\
              "functional.mean_conc   = {0.mean_conc}\n"\
              "functional.stddev_conc = {0.stddev_conc}\n"\
              "functional.cm_rate     = {0.cm_rate}\n"\
              .format(self, repr(self.rootworkdir), repr(self.age))
    return header + string
