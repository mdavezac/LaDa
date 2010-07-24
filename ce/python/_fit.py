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
    from lada import ce 
    import numpy


    self._pis       = None
    self._energies  = None
    self._weights   = None
    self._str_onoff = None
    self._cls_onoff = numpy.array([True for u in clusters], dtype=bool)
    self._classes = ce.MLClusterClasses(clusters)
    self._structures = None

  def __call__(self, A = None, b = None):
    """ Returns observation matrix A and target vector b.

        if A or b are given, adds to the pre-existing matrix and vector.  A and
        b can then be used in linear-least square fit method
        numpy.linalg.lstsq(A,b)
    """

    import numpy
    import pyublas

    # gets array lengths
    x, y = 0, 0
    for u in self._str_onoff:
      if u : x += 1
    for u in self._cls_onoff:
      if u : y += 1

      
    # creates arrays if they are not provided.
    if A == None: A = numpy.zeros((x,y), dtype='float64')
    if b == None: b = numpy.zeros((x,), dtype='float64')

    # checks arrays have correct size
    assert A.shape == (x, y), "A does not have correct shape.\n"
    assert len(b) == x, "b does not have correct length.\n" 

    # fills arrays from local data.
    mat = self._pis[self._str_onoff,:][:,self._cls_onoff].copy()
    vec = self._energies[self._str_onoff].copy()
    weights = self._weights[self._str_onoff]
    for i in range(len(b)): mat[i,:] *= weights[i] 
    vec *= weights

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
    if self._pis == None:  # on-offs for clusters.
      nclasses = len(self._classes)
      self._pis = numpy.zeros((1, nclasses), dtype='float64').copy()
      self._energies = numpy.zeros((1,), dtype='float64').copy()
      self._weights = numpy.zeros((1,), dtype='float64').copy()
      self._str_onoff = numpy.array([True], dtype=bool)
      self._structures = [structure]
    else: 
      self._energies.resize((len(self._energies)+1,))
      self._weights.resize((len(self._weights)+1,))
      self._str_onoff.resize((len(self._str_onoff)+1,))
      self._pis.resize((self._pis.shape[0]+1, self._pis.shape[1]))
      self._structures.append(structure)

    # computes pis and adds data
    self._pis[-1,:] = self._classes.pis(structure)
    self._energies[-1] = structure.energy
    self._weights[-1] = structure.weight
    self._str_onoff[-1] = True

  def extinguish_structure(self, n = None, on = "all"):
    """ Extinguishes structure at index n (n can be a list of indices).
        L{Structure} at index 0 is the first structure added to the set.
        L{Structure} at index 1 is the second structure added to the set.
        And so on...
        If on == "all", then turns all other structures on.
        If n is None and on == all, then  all structures are turned on.
    """

    if on == "all": self._str_onoff[:] = True
    if n != None: self._str_onoff[n] = False

  def extinguish_cluster(self, n = None, on = "all", mask=False):
    """ Extinguishes cluster class at index n (n can be a list of indices).
        If on == "all", then turns all other cluster classes on.
        If n is None and on == all, then  all cluster classes are turned on.
    """

    if on == "all": self._cls_onoff[:] = not mask
    if n != None: self._cls_onoff[n] = mask

  def assign_genome(self, n):
    """ Assigns a bitstring genome for fitting.
        A bitstring genome should be a numpy boolean array of the same size as
        the number of cluster classes.
    """

    if hasattr(n, "shape"):
      assert len(n.shape) == 1, "Wrong shape for parameter n.\n" 
    self._cls_onoff[:] = n

  def set_weight(self, n, weight): 
    """ Sets the weigh of structure at index n (n can be a list of indices.) """
    if hasattr(n, "shape"):
      assert len(n.shape) == 1, "Wrong shape for parameter n.\n" 
    self._weights[n] = weight

  def size(self): 
    """ Returns the total number of cluster and structures (whether on or off). """
    return len(self._cls_onoff), len(self._str_onoff)


  def reset_clusterclasses(self, classes):
    """ Resets cluster classes to those on input.
        This is necessary after changing the cluster classes. 
        The weights and energies remain the same after and before call.
    """
    from numpy import array
    from lada.ce import MLClusterClasses

    self._cls_onoff = array([True for u in classes], dtype=bool)
    self._classes = MLClusterClasses(classes)

    structures = self._structures
    weights = self._weights
    self._pis = None
    for structure in structures:
      self.add_structure(structure)
    self._weights = weights


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
    


class PairRegulatedFit(Fit):
  """ CE with pair regularisation.
      This class returns (__call__) observation and target arrays for use with
      linear least square fit routine from numpy.linalg.
      The last observations are for regularisation.
      For some reason (a bug that lasted long enough to become a feature?)
      there are two types of reugalarization used at NREL. They differ only by
      in the exact normalization (or equivalently, they differ in the unit of tcoef).
      self.alpha and self.tcoef parameterize the regularization.
  """

  def __init__(self, clusters, alpha = 2, tcoef = 1, type="laks"):
    """ alpha and tcoef are the regularization coefficients. """
    from numpy.linalg.basic import norm as np_norm

    # constructs from base.
    Fit.__init__(self, clusters)

    # now adds regulatization.
    self.type = type
    self.alpha = alpha
    self.tcoef = tcoef

    # and initializes regularization data.
    self._norms = []
    for i, class_ in enumerate(self._classes): 
      if class_.order != 2: continue
      self._norms.append( (i,np_norm(class_[0].origin.pos - class_[0][0].pos)**2) )
      
  def reset_clusterclasses(self, classes):
    """ Recompute pis from known structures.

        This is necessary after changing the cluster classes. 
        The weights and energies remain the same after and before call.
        Bookkeeping for regulated pairs is re-initialized.
    """
    from numpy.linalg.basic import norm as np_norm
    Fit.reset_clusterclasses(self, classes)
    self._norms = []
    for i, class_ in enumerate(self._classes): 
      if class_.order != 2: continue
      self._norms.append( (i,np_norm( class_[0].origin.pos - class_[0][0].pos)) )
      

  def _compute_weights(self, A):
    from math import pow, sqrt, fabs
     
    # Get data for "on" pairs.
    pairs = []
    for p, norm in self._norms:
      if self._cls_onoff[p]: pairs.append((p,norm))

    assert len(pairs) == A.shape[0], "Wrong size for A."

    weights = []
    normalization = 0e0

    # Constructs coefficients
    coef = 4e0
    if self.type != "laks": coef = 1e0
    for i, norm in pairs: 
      w = pow( norm, self.alpha ) * 0.5e0
      weights.append( (i, sqrt(w)) )
      normalization += w * coef

    if self.type != "laks": normalization = sqrt(fabs(self.tcoef)) / normalization
    else: normalization = sqrt(fabs(self.tcoef) / normalization)

    # now constructs matrix.
    for i, data in enumerate(weights):
      j = 0
      for u in range(data[0]):
        if self._cls_onoff[u]: j+=1
      A[i, j] = data[1] * normalization

  def __call__(self, A = None, b = None):
    """ Returns observation matrix A and target vector b.
        if A or b are given, adds to the pre-existing matrix and vector.
        A and b can then be used in linear-least square fit method numpy.linalg.lstsq(A,b)
        A and b (if given) should be large enough to accomodate for pair regularization.
    """
    from numpy import zeros
    from math import fabs

    if fabs(self.tcoef) < 1e-12: return Fit.__call__(self, A, b)

    # gets array lengths
    x, y = 0, 0
    for u in self._str_onoff:
      if u : x += 1
    for u in self._cls_onoff:
      if u : y += 1
    xreg = x
    for p, norm in self._norms:
      if self._cls_onoff[p]: xreg += 1

      
    # creates arrays if they are not provided.
    if A == None: A = zeros((xreg,y), dtype='float64')
    if b == None: b = zeros((xreg,), dtype='float64')


    # calls base class
    if xreg == x:
      return Fit.__call__(self, A[:x,:], b[:x]) 
    A[:x,:], b[:x] = Fit.__call__(self, A[:x,:], b[:x]) 

    # adds regularization
    self._compute_weights(A[x:,:])

    return A, b 

def leave_one_out( fitter ):
  """ Performs leave-many out on fitter using the given sets. 
      Returns a matrix where each line contains the fitted and predicted structures.
      The predictions are the diagonal elements.
  """
  import numpy

  ncls, nstr = fitter.size()  

  # matrix which holds errors for all structures. 
  errors = numpy.zeros( (nstr,nstr), dtype='float64')
  # fitting matrices from which to get errors.
  fitter.extinguish_structure() 
  A_all, b_all = fitter()
  A_all = A_all[:nstr, :].copy()
  b_all = b_all[:nstr].copy()

  # loop over fitting sets.
  for i in xrange(nstr):
    fitter.extinguish_structure(i)
    A, b = fitter()
    x, residues, rank, s = numpy.linalg.lstsq(A, b)
    errors[i, :] = numpy.dot(A_all, x) - b_all

  return errors


def leave_many_out( fitter, sets ):
  """ Performs leave-many out on fitter using the given sets. 
      Returns a matrix where each line contains the fitted and predicted structures for one set.
  """
  import numpy

  ncls, nstr = fitter.size()  

  # matrix which holds errors for all structures. 
  errors = numpy.zeros( (len(sets),nstr), dtype='float64')
  # fitting matrices from which to get errors.
  fitter.extinguish_structure() 
  A_all, b_all = fitter()
  A_all = A_all[:nstr, :].copy()
  b_all = b_all[:nstr].copy()

  # loop over fitting sets.
  for i, s in enumerate(sets):
    fitter.extinguish_structure(s)
    A, b = fitter()
    x, residues, rank, s = numpy.linalg.lstsq(A, b)
    errors[i, :] = numpy.dot(A_all, x) - b_all

  return errors

