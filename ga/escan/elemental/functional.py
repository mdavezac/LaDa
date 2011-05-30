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
        :Kwarg nmax: 
          Maximum size of bitstrings. Defaults to 20.
        :Kwarg nmin: 
          Minimum size of bitstrings.
        :Kwarg crossover_rate:
          Rate for crossover operation
        :Kwarg swap_rate:
          Rate for specie-swap mutation-operation
        :Kwarg growth_rate:
          Rate for the mutation operator where layers are inserted or removed.
        :Kwarg rate: 
          Offspring rate. Defaults to ``max(0.1, 1/popsize)``.
        :Kwarg popsize: 
          Population size. Defaults to 100.
        :Kwarg maxgen: 
          Maximum number of generations. Defaults to 0.
        :Kwarg current_gen:
          Current generation. Defaults to 0 or to whatever is in the restart.
    """
    from .. import CompareSuperCells
    from .operators import Crossover, GrowthMutation, SwapMutation
    from ...standard import Mating
    super(Darwin, self).__init__(evaluator, **kwargs)

    # add parameters from this GA.
    self.nmin = kwargs.pop('nmin', -1)
    """ Minimum bistring size. """
    self.nmax = kwargs.pop('nmax', 20)
    """ Maximum bistring size. """

    # creates mating operation.
    self.crossover_rate = kwargs.pop('crossover_rate', 0.8)
    """ Rate for crossover over other operations. """
    self.swap_rate = kwargs.pop('swap_rate', 0.8)
    """ Rate for swap-like mutations over other operations. """
    self.growth_rate = kwargs.pop('growth_rate', 0.8)
    """ Rate for growth-like mutations over other operations. """
    self.matingops = Mating(sequential=False)
    if self.crossover_rate > 0e0: self.matingops.add(Crossover(self.nmin, self.nmax, step=2))
    if self.swap_rate > 0e0:      self.matingops.add(SwapMutation(-1))
    if self.growth_rate > 0e0:    self.matingops.add(GrowthMutation(self.nmin, self.nmax, step=2))

  def mating(self):
    """ Calls mating operations. """
    return self.matingops(self)

  def compare(self, a, b):
    """ Compares two bitstrings. """
    from numpy import all
    if len(a.genes) != len(b.genes): return False
    return all( a.genes == b.genes )

  def taboo(self, indiv):
    """ No two individuals in the population and the offspring are the same. """
    from itertools import chain
    assert len(indiv.genes) >= self.nmin and len(indiv.genes) <= self.nmax,\
           (indiv.genes, len(indiv.genes), self.nmin, self.nmax)
    for a in chain(self.population, self.offspring):
      if self.compare(a, indiv): return True
    return False

  def Individual(self):
    """ Generates an individual (with mpi) """
    from . import Individual
    result = Individual(nmax=self.nmax, nmin=self.nmin)
    assert len(result.genes) >= self.nmin and len(result.genes) <= self.nmax
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
              .format(self, repr(self.rootworkdir), repr(self.age))
    return header + string
