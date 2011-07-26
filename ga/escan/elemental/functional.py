""" GA as a functional, for use with ladajobs. """
__docformat__ = "restructuredtext en"
from ...functional import Darwin as DarwinBase

__all__ = ['Darwin']

class Darwin(DarwinBase): 
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
        :Kwarg dosyum:
          Whether to look for inversion/translation symmetries using the simplistic
          implementation of Individual.
    """
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
    """ Container/Functor which calls the mating operations. """
    if self.crossover_rate > 0e0:
      self.matingops.add(Crossover(self.nmin, self.nmax, step=2), self.crossover_rate)
    if self.swap_rate > 0e0:     
      self.matingops.add(SwapMutation(-1), self.swap_rate)
    if self.growth_rate > 0e0:   
      self.matingops.add(GrowthMutation(self.nmin, self.nmax, step=2), self.growth_rate)
    if 'dosym' in kwargs:
      self.dosym = kwargs.pop('dosym')
      """ Whether or not to use simplistic symmetries when comparing individuals. """

  def mating(self):
    """ Calls mating operations. """
    return self.matingops(self)


  def Individual(self):
    """ Generates an individual (with mpi) """
    from . import Individual
    result = Individual(nmax=self.nmax, nmin=self.nmin, dosym=getattr(self, 'dosym', False))
    assert len(result.genes) >= self.nmin and len(result.genes) <= self.nmax
    return result

  def __repr__(self):
    """ Returns representation of this instance. """
    header = "from {0.__class__.__module__} import {0.__class__.__name__}\n"\
             "from {0.comparison.__module__} import {0.comparison.__name__}\n"\
             "from {0.history.__module__} import {0.history.__name__}\n".format(self)
    string = repr(self.evaluator) + "\n" + repr(self.matingops) + "\n"
    string += "functional = {0.__class__.__name__}(evaluator)\n"\
              "functional.pools       = {0.pools}\n"\
              "functional.nmax        = {0.nmax}\n"\
              "functional.nmin        = {0.nmin}\n"\
              "functional.popsize     = {0.popsize}\n"\
              "functional.max_gen     = {0.max_gen}\n"\
              "functional.current_gen = {0.current_gen}\n"\
              "functional.rate        = {0.rate}\n"\
              "functional.current_gen = {0.current_gen}\n"\
              "functional.comparison  = {0.comparison.__name__}\n"\
              "functional.age         = {1}\n"\
              "functional.history     = {2}\n"\
              .format(self, repr(self.age), repr(self.history))
    if getattr(self, 'dosym', False): 
      string += "functional.dosym       = {0.dosym}\n".format(self)
    string += "# functional.rootworkdir, please set evaluator.outdir instead\n"
    return header + string
