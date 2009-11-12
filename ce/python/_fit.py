#
#  Version: $Id$
#


class Fit(): 
  """ Computes lstsq vector and matrix for a sub-set of clusters and a sub-set of structures.
  """
  def __init__(self, clusters):
    """ Initialises the fitting class.
        clusters should be a ce.MLClusterClass instance, or convertible. A copy
        of clusters is assigned to self.classes
    """
    import copy
    from lada import ce 


    self._data = None
    self._classes = clusters

  def __call__(self, A = None, b = None):
    """ Returns observation matrix A and target vector b.
        if A or b are given, adds to the pre-existing matrix and vector.
        A and b can then be used in linear-least square fit method numpy.linalg.lstsq(A,b)
    """

    import numpy.linalg
    import pyublas

    x, y = 0, 0
    for u in self._data[1:,-1]:
      if u : x += 1
    for u in self._data[0,:-3]:
      if u : y += 1

    assert A == None or A.shape == (x, y), "A does not have correct shape.\n"
    assert b == None or len(b) == x, "b does not have correct length.\n" 
      
    if A == None and b == None:
      A = numpy.narray([0e0 for i in range(x*y)]).shape(x,y)
      b = numpy.narray([0e0 for i in range(y)]).shape(1,y)
    elif A == None: 
      x, y = 0, len(b) 
      A = numpy.narray([0e0 for i in range(x*y)]).shape(x,y)
    elif A == None: 
      y = A.shape[0]
      b = numpy.narray([0e0 for i in range(y)]).shape(1,y)

    str_ons = self._data[:,-1]
    cls_ons = self._data[0,:]

    weights = self._data[str_on,-2]
    mat = self._data[str_ons, cls_ons].copy()
    vec = self._data[str_ons, -3].copy()
    for i in range(len(b)):
      mat[i] *= w[i] 
      vec[i,:] *= w[i]

    return A+mat, b+vec



  def add_structure(self, structure):
    """ Adds a structure to input set.
        Computes the pis for given set of clusters.
        Expects structure.energy to hold the target energy, 
        and structure.weight the weight of the fitting set in the structure.
    """
    from lada.ce import find_pis
    import numpy
    import pyublas
    
    # resize array to fit structure.
    if self._data == None:  # on-offs for clusters.
      nclasses = len(self._classes)
      self._data = numpy.array( [float(1) for o in xrange(nclasses+3)] ).reshape(1, nclasses+3)
      self._data[:] = True
      self._data[-3:] = False

    self._data = numpy.resize(self._data, (self._data.shape[0]+1, self._data.shape[1]))
    print self._data

    # computes pis and adds data
    self._data[-1,:-3] = self._classes.pis(structure)
    self._data[-1,-3] = structure.energy
    self._data[-1,-2] = structure.weight
    self._data[-1,-1] = True

  def extinguish_structure(self, n = None, on = "all"):
    """ Extinguishes structure at index n (n can be a list of indices).
        Structure at index 0 is the first structure added to the set.
        Structure at index 1 is the second structure added to the set.
        And so on...
        If on == "all", then turns all other structures on.
        If n is None and on == all, then  all structures are turned on.
    """

    if on == "all": self._data[1:, -1] = True
    if n != None: self._data[n+1, -1] = False

  def extinguish_cluster(self, n = None, on = "all"):
    """ Extinguishes cluster class at index n (n can be a list of indices).
        If on == "all", then turns all other cluster classes on.
        If n is None and on == all, then  all cluster classes are turned on.
    """

    if on == "all": self._data[0, :-3] = True
    if n != None: self._data[0, n] = False
 
  def assign_genome(self, n):
    """ Assigns a bitstring genome for fitting.
        A bitstring genome should be a numpy boolean array of the same size as
        the number of cluster classes.
    """

    if hasattr(n, "shape"):
      assert len(n.shape) == 1, "Wrong shape for parameter n.\n" 
    if hasattr(n, "__len__"):
      assert len(n) == self._data.shape[0]-3, "n too long.\n" 
    self._data[0, :-3] = n

  def set_weight(self, n, weight): 
    """ Sets the weigh of structure at index n (n can be a list of indices.) """
    self._data[n, -2] = weight

  def size(self): 
    """ Returns the total number of cluster and structures (whether on or off). """
    return self._data.shape[0]-1, self._data.shape[1]-3



  def read_directory(self, path):
    """ Reads a directory containing LDAs.dat file and adds the structures to the fitting set. 
        The structure files should exist in the same directory.
    """
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
    


     

    
    


