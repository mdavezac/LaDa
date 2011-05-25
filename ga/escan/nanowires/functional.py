""" GA functional for Nanowire optimization. """
__docformat__ = "restructuredtext en"
__all__ = ['Darwin']

# from ..elemental.extract import Extract as GAExtract
from ..functional import Darwin as ElementalDarwin

class Darwin(ElementalDarwin):
  """ GA functional for optimizations of epitaxial structures. """

  def __init__(self, evaluator, **kwargs): 
    """ Initializes a GA functional. 
         
        :Parameters:
          evaluator : `lada.ga.escan.elemental.Bandgap`
            Functional which uses vff and escan for evaluating physical properties.
        :Kwarg Nmax: 
          Maximum size of bitstrings. Defaults to 20.
        :Kwarg Nmin: 
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
    from ...standard import Mating
    from ...bitstring import VariableSizeCrossover, SwapMutation, GrowthMutation
    # add parameters from this GA.
    self.Nmin = kwargs.pop('Nmin', -1)
    """ Minimum bistring size. """
    self.Nmax = kwargs.pop('Nmax', 20)
    """ Maximum bistring size. """

    # calls parent constructor.
    super(Darwin, self).__init__(evaluator, **kwargs)

    # creates mating operation.
    self.crossover_rate = kwargs.pop('crossover_rate', 0.8)
    """ Rate for crossover over other operations. """
    self.swap_rate = kwargs.pop('swap_rate', 0.8)
    """ Rate for swap-like mutations over other operations. """
    self.growth_rate = kwargs.pop('growth_rate', 0.8)
    """ Rate for growth-like mutations over other operations. """
    self.mating = Mating(sequential=False)
    if self.crossover_rate > 0e0: self.mating.add(VariableSizeCrossover(self.Nmin, self.Nmax))    
    if self.swap_rate > 0e0:      self.mating.add(SwapMutation(-1))
    if self.growth_rate > 0e0:    self.mating.add(GrowthMutation(self.Nmin, self.Nmax))    

  def compare(self, a, b):
    """ Compares two bitstrings. """
    from numpy import all
    if len(a.genes) != len(b.genes): return False
    return all( a.genes == b.genes )

  def Individual(self):
    """ Generates an individual (with mpi) """
    from . import Individual
    return Individual(Nmax=self.Nmax, Nmin=self.Nmin)

  def __repr__(self):
    """ Returns representation of this instance. """
    raise NotImplementedError("Cannot represent this functional.")
