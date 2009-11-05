#
#  Version: $Id$
#


class Fit(): 
  def __init__(self, clusters):
    """
      Clusters should be an array with the following form:
        [ centered on Sublattice 0, centered on Sublattice 1, ... ]
    """
    import copy


    self.clusters = copy.deepcopy(clusters)
    self.pis = numpy.array([])
    self.energies = numpy.array([])
    self.structures = []
    self.weights = []

  def add_structure(structure):
    from lada.ce import find_pis
    
    


