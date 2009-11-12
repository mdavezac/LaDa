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

  def __call__(self, A = None, b = None, rcond = -1): 

    import numpy.linalg
    import pyublas

    x, y = 0, 0
    for u in self._data.shape[1:,-1]): if u > 0e0: x += 1
    for u in self._data.shape[0,::-1]): if u > 0e0: y += 1

    assert A == None or A.shape == x, y, "A does not have correct shape.\n"
    assert b == None or len(b) == x, "b does not have correct length.\n" 
      
    if A == None and b == None:
      A = numpy.narray([0e0 for i in range(x*y)]).shape(x,y)
      b = numpy.narray([0e0 for i in range(y)]).shape(1,y)
    elif A == None: 
      x, y = 0, len(b) 
      for u in self._data.shape[1:,-1]): if u > 0e0: x += 1
      A = numpy.narray([0e0 for i in range(x*y)]).shape(x,y)
    elif A == None: 
      y = A.shape[0]
      b = numpy.narray([0e0 for i in range(y)]).shape(1,y)

    str_ons = self._data[:,-1]
    cls_ons = self._data[0,:]

    weights = self._data[str_on,-2]
    A = self._data[str_ons, cls_ons].copy()
    b = self._data[str_ons, -3].copy()
    for i in range(len(b)):
      b[i] *= w[i] 
      A[i,:] *= w[i]

    return A, b



  def add_structure(self, structure):
    from lada.ce import find_pis
    import numpy
    import pyublas
    
    nclasses = len(self.classes)
    # resize array to fit structure.
    if self._data == None:  # on-offs for clusters.
      self._data = numpy.array( [float(1) for o in xrange(nclasses+3)] ).reshape(1, nclasses+1)
      for u in range(len(self._data)): self._data[u] = True
      self._data[-3:] = False

    self._data = numpy.resize(self._data, (self._data.shape[0]+1, self._data.shape[1]))

    # computes pis and adds data
    self._data[-1,:nclasses] = self.classes.pis(structure)
    self._data[-1,-3] = structure.energy
    self._data[-1,-2] = structure.weight
    self._data[-1,-1] = True

  def extinguish_structure(self, n, on = "all"):

    if on == "all": self._data[1:, -1] = True
    self._data[n+1, -1] = False

  def extinguish_cluster(self, n, on = "all"):

    if on == "all": self._data[0, :-3] = True
    self._data[0, n] = False
 
  def assign_genome(self, n, on = "all"):

    if hasattr(n, "shape"):
      assert len(n.shape) == 1, "Wrong shape for parameter n.\n" 
    if hasattr(n, "__len__"):
      assert len(n) == self._data.shape[0]-3, "n too long.\n" 
    self._data[0, :-3] = n

  def set_weight(self, n, weight): self._data[n, -2] = weight

  def size(self): return self._data.shape[0]-1, self._data.shape[1]-3



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
    print self._data 
    


     

    
    


