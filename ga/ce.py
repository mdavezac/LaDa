#
#  Version: $Id: ce.py 1396 2009-11-16 02:59:48Z davezac $
#

class Eval:

  def __init__( self, lattice, path="data", alpha=2, tcoef=10, lmo_ratio = 0.333333, **keywords):
    """ Creates an evaluation function for a multi-site cluster expansion.
        Requires on input the lattice, a path towards a directory of structure data,
        and the a short form for the allowed clusters.
    """
    from lada import ce 
    import os
    import re
    import random

    self.lattice = lattice

    # data for the fit.
    self.alpha, self.tcoef = alpha, tcoef

    # checks that keywords are well formed.
    key_regex = re.compile("B(\d+)")
    for key in keywords.keys(): 
      a_re = re.match(key_regex, key)
      assert a_re != None, "Keyword %s is not of the form B(\d+)" % (key)
      assert a_re.end() == len(key), "Keyword %s is not of the form B(\d+)" % (key)
      assert int(a_re.group(1)) > 1, "Cannot understand keyword %s" % (key)
      assert keywords[key] > 0, "Cannot understand input %s=%i" % (key, keywords[key])

    # sanity check.
    assert len(keywords) > 0, "No keywords found on input. Can't construct clusters.\n"

    keys = sorted(keywords.keys(), cmp)

    def equivalent_sites():
      """ Returns a list containing only one site index per inequivalent sub-lattice. """
      from lada import crystal
    
      result = set( i for i in range(len(self.lattice.sites)) ) 
      for i, site in enumerate(self.lattice.sites):
        if i not in result: continue
        for op in lattice.space_group:
          j = crystal.which_site( op(site.pos), self.lattice )
          if j != i and j in result: result.remove(j)
      return result

    # now creates multi-lattice clusters, along with index bookkeeping.
    self._mlclasses = ce.create_clusters(lattice, nth_shell=0, order=0, site=0) # J0
    # creates J1 for these sites.
    for site in equivalent_sites():
      self._mlclasses.extend(ce.create_clusters(lattice, nth_shell=0, order=1, site=site))
    # creates many bodies.
    self._fixed = len(self._mlclasses)
    for site in equivalent_sites():
      for key in keys:
        regex = re.match(key_regex, key)
        order = int(regex.group(1))
        shell = keywords[key]
        self._mlclasses.extend(ce.create_clusters(lattice, nth_shell=shell, order=order, site=site))

    # creates fitting function 
    self.fitter = ce.PairRegulatedFit(self._mlclasses, alpha=self.alpha, tcoef=self.tcoef)
    self.fitter.read_directory(path)

    # now reads/creates set.
    set_path = os.path.join(path, "set")
    self._sets = []
    ncls, nstr = self.fitter.size()
    if os.path.exists(set_path):
      # reads path 
      print "Found input set %s:\n" % (set_path)
      file = open(set_path, "r")
      for line in file:
        self._sets.append( [int(i) for i in line.spli()] )
        print "  ", self._sets[-1]
        for i in self._sets[-1]:
          assert i < nstr, "Index in set is larger than the number of structures.\n"
      print
    else: 
      # creates set.
      size = max( int(float(nstr)*lmo_ratio), 1 )
      if size == 1: "Will perform leave-one out cross-validation.\n"
      else: 
        print "Creating leave-many out set:\n"
        full_set = []
        for i in xrange(size):
          if len(full_set) < size:
            full_set = [ j for j in range(nstr) ]
            random.shuffle(full_set)
          self._sets.append( full_set[:size] )
          full_set = full_set[size:]
          print "  ", self._sets[-1]

  def __len__(self):
    """ Returns the number of many-body figures (minus J0, J1, and pairs) """
    return len(self._mlclasses) - self._fixed
     


  def __call__(self, indiv):
    """ Returns cv score for this individual.
        Optimizes cv score with respect to pairs, and pair regularization.
    """
    from math import fabs
    from .. import ce
    import numpy

    assert len(self._mlclasses) - self._fixed == len(indiv.genes),\
           "Individual is of incorrect size.\n"
    onoffs = numpy.zeros((len(self._mlclasses),), dtype=bool)
    onoffs[:self._fixed] = True
    onoffs[self._fixed:] = indiv.genes

    self.fitter.extinguish_cluster(onoffs, mask=True) 
    self.fitter.alpha = self.alpha
    self.fitter.tcoef = self.tcoef

    if len(self._sets) == 0:
      errors = ce.leave_one_out(self.fitter)
      npreds, nstr = errors.shape
      prediction = errors[ numpy.arange(errors.shape[0]), numpy.arange(errors.shape[0]) ]
      return numpy.average(prediction*prediction)
    else:
      errors = ce.leave_many_out(self.fitter, self._sets)
      npreds, nstr = errors.shape
      prediction = errors[ [[(j,i) for i in self._sets[j]] for j in len(self.sets)] ]
      return numpy.average(prediction*prediction)

class Individual:
  """ An individual for boolean bitstrings.
      The size of the individual is given statically by Individual.size
  """

  size = 10

  def __init__(self):
    """ Initializes a boolean bitstring individual randomly.
    """
    from random import choice
    import numpy

    self.genes = numpy.array([ choice([True,False]) for i in xrange(Individual.size) ], dtype=bool)
    self.evaluated = False
  
  def __eq__(self, a): 
    from math import fabs
    if a == None: return False
    for i, b in enumerate(a.genes):
      if self.genes[i] != b: return False
    return True
