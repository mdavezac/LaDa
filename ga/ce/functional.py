""" GA as a functional, for use with ladajobs. """
__docformat__ = "restructuredtext en"
__all__ = ['Darwin']
from ..functional import Darwin as StandardDarwin

class Darwin(StandardDarwin):
  """ Main GA functional for CE. """
  def __init__(self, evaluator, **kwargs):
    """ Initializes a GA functional. 
         
        :Parameters:
          evaluator : `lada.ga.escan.elemental.Bandgap`
            Functional which uses vff and escan for evaluating physical properties.
        :Kwarg rate: 
          Offspring rate. Defaults to ``max(0.1, 1/popsize)``.
        :Kwarg popsize: 
          Population size. Defaults to 100.
        :Kwarg maxgen: 
          Maximum number of generations. Defaults to 0.
        :Kwarg current_gen:
          Current generation. Defaults to 0 or to whatever is in the restart.
        :Kwarg pools:
          Number of pools over which to perform calculations. 
        :Kwarg rootworkdir:
          Root of the working directory where actual calculations are carried
          out. Defaults to $SCRACTH
    """
    super(Darwin, self).__init__(evaluator, **kwargs)
    
