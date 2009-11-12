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
    self._pis = None
    self._energies = None
    self._onofss = None
    self.structures = []

  def add_structure(self, structure):
    from lada.ce import find_pis
    import numpy
    import pyublas
    
    nclasses = len(self.classes)
    # resize array to fit structure.
    if self._data == None: 
      self._data = numpy.array( [float(0) for o in xrange(nclasses+1)] ).reshape(1, nclasses+1)
    else:
      self._data = numpy.resize(self._data, (self._data.shape[0]+1, self._data.shape[1]))

    # computes pis and assigns them.
    self._data[-1,:len(self.classes)] = self.classes.pis(structure)
    # records energy.
    self._data[-1,-1] = structure.energy
    print self._data[-1]

  def read_directory(self, path):
    import os.path
    import re
    from lada import crystal

    assert os.path.exists(path), "%s does not exist.\n" % (path)
    assert os.path.isdir(path), "%s is not a directory.\n" % (path)
    LDAs_filename = os.path.join(path, "LDAs.dat")
    assert os.path.exists(LDAs_filename), "%s does not exist.\n" % (LDAs_filename)

    LDAs = open(LDAs_filename, "r")
    for line in LDAs:
      if re.match("\s*#", line) != None: continue 
      filename = line.split()[0] 
      structure = crystal.read_structure( os.path.join(path, filename) )
      structure.energy = float(line.split()[1])
      structure.weight = 1e0
      self.add_structure(structure)
    


     

    
    


