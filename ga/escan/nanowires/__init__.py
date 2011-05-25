""" A GA subpackage defining standard genetic operator for nanowires. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Individual']

from ...bitstring import Individual as BitstringIndividual
class Individual(BitstringIndividual):
  """ An individual for nanowires.
      
      Comparison between two individuals expect that a bitstring represents an
      elemental nanowires.
  """
  def __init__(self, Nmin=0, Nmax=20):
    """ Initializes a bitstring individual randomly. """
    from random import randint
    import numpy

    self.size = randint(Nmin, Nmax)
    super(Individual, self).__init__(self.size)
    self.genes = numpy.array([ int(randint(0,1)) for i in xrange(self.size) ], dtype="int")
  
  def __eq__(self, a): 
    """ Compares two elemental nanowires. """
    from numpy import all
    if a == None: return False
    if len(a.genes) != len(self.genes): return False
    return all(self.genes == a.genes)

def exec_input(script, namespace = None):
  """ Executes an input script including namespace for ga nanowires. """ 
  from ....escan import exec_input as escan_exec_input
  from functional import Darwin
  from converter import Converter
  from ..elemental.evaluator import Bandgap as BandGapEvaluator, Dipole as DipoleEvaluator

  dictionary = {}
  if namespace != None: dictionary.update(namespace)
  dictionary['Individual'] = Individual
  dictionary['Darwin'] = Darwin
  dictionary['Converter'] = Converter
  dictionary['BandGapEvaluator'] = BandGapEvaluator
  dictionary['DipoleEvaluator']  = DipoleEvaluator
  return escan_exec_input(script, dictionary)

def read_input(filepath = "input.py", namespace = None):
  """ Reads an input file including namespace for ga nanowires. """ 
  from ....escan import read_input as escan_read_input
  from functional import Darwin
  from converter import Converter
  from ..elemental.evaluator import Bandgap as BandGapEvaluator, Dipole as DipoleEvaluator

  dictionary = {}
  if namespace != None: dictionary.update(namespace)
  dictionary['Individual'] = Individual
  dictionary['Darwin'] = Darwin
  dictionary['Converter'] = Converter
  dictionary['BandGapEvaluator'] = BandGapEvaluator
  dictionary['DipoleEvaluator']  = DipoleEvaluator
  return escan_read_input(filepath, dictionary)
