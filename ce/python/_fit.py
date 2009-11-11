#
#  Version: $Id$
#


class Fit(): 
  """ Fits a sub-set of clusters to a sub-set of structures.
      _ self.classes is the full set of clusters. It is an instance of ce.MLClusterClasses.
      _ self.onoffs is a vector specifying which clusters are to be used.
      _ self.structures is a set of structures.
      
  """
  def __init__(self, clusters):
    """ Initialises the fitting class.
        clusters should be a ce.MLClusterClass instance, or convertible. A copy
        of clusters is assigned to self.classes
    """
    import copy
    from lada import ce 


    self.classes = ce.MLClusterClasses(clusters)
    self._data = None
    self.structures = []
    self.onoffs = [] 

  def add_structure(structure):
    from lada.ce import find_pis
    import numpy
    import pyublas
    
    nclasses = len(self.classes)
    # resize array to fit structure.
    if self._data == None: 
      self._data = numpy.arange(nclasses+1).reshape(1, nclasses+1)
    else: numpy.resize(self._data, (self.shape[0]+1, self.shape[1]))

    # computes pis and assigns them.
    self._data[-1,:nself.classes] = ce.find_pis(self.classes, structure)
    # records energy.
    self._data[-1,-1] = structure.energy

     

    
    


