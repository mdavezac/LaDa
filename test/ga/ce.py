#
#  Version: $Id$
#

class Eval:

  def __init__( self, lattice, path="data",
                alpha=(0,0), tcoef=(0,0), pairs=(1,20),
                lmo_ratio = 0.333333, **keywords):
    """ Creates an evaluation function for a multi-site cluster expansion.
        Requires on input the lattice, a path towards a directory of structure data,
        and the a short form for the allowed clusters.
    """
    from lada import ce 
    import re

    self.lattice = lattice

    # data for the fit.
    assert len(pairs) == 2, "value of \"pairs\" keyword should contain min and max allowed pairs."
    self.pairs = min(pairs), max(pairs)
    assert pairs[0] > 0 and pairs > 0, "Negative value for pairs is not allowed."
    assert len(alpha) == 2, "value of \"alpha\" keyword should contain min and max allowed regularization alpha."
    assert alpha[0] >= 0e0 and alpha[1] >= 0e0, "Negative value for alpha is not allowed."
    self.alpha = min(alpha), max(alpha)
    assert len(tcoef) == 2, "value of \"tcoef\" keyword should contain min and max allowed regularization tcoef."
    assert tcoef[0] >= 0e0 and tcoef[1] >= 0e0, "Negative value for tcoef is not allowed."
    self.tcoef = min(tcoef), max(tcoef)
    assert "alpha" not in keywords.keys(), "WTF"

    # checks that keywords are well formed.
    for key in keywords.keys(): 
      a_re = re.match("(\d+)B", key)
      assert a_re != None, "Keyword %s is not of the form J[0,1] or (\d+)B(\d+)" % (a)
      assert a_re.end() == len(a), "Keyword %s is not of the form J[0,1] or (\d+)B(\d+)" % (a)
      assert int(a_re.group(1)) > 2, "Cannot understand keyword %s" % (a)
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

    # now creates multi-lattice clusters.
    self.mlclasses = ce.create_clusters(lattice, nth_shell=0, order=0, site=site) # J0
    for site in equivalent_sites():
      # creates pairs.
      self.mlclasses.extend(ce.create_clusters(lattice, nth_shell=self.pairs[1], order=2, site=site))
      # creates J1 for these sites.
      self.mlclasses.extend(ce.create_clusters(lattice, nth_shell=0, order=1, site=site))
      # creates many bodies.
      for key in keys:
        regex = re.compile("(\d+)B")
        order = int(regex.group(1))
        shell = keywords[key]
        self.mlclasses.extend(ce.create_clusters(lattice, nth_shell=shell, order=order, site=site))

    # creates fitting function 
    self.fitter = ce.RegulatedFit(mlclasses, alpha=self.alpha[0], tcoef=self.tcoef[0])
    self.fitter.read_directory(path)

    # now reads/creates set.
    set_path = os.path.join(directory, "set")
    self.sets = []
    ncls, nstr = self.fitter.size()
    if os.path.exists(set_path):
      # reads path 
      print "Found input set %s:\n" % (set_path)
      file = open(set_path, "r")
      for line in file:
        self.sets.append( [int(i) for i in line.spli()] )
        print "  ", self.sets[-1]
        for in self.sets[-1]:
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
          self.sets.append( full_set[:size] )
          full_set = full_set[size:]
          print "  ", self.sets[-1]

  def __len__(self):
    """ Returns the number of many-body figures (minus J0, J1, and pairs) """
    return len(self.mlclasses)
     
  def _fit(self):

    if len(self.sets) == 0:
      errors = ce.leave_one_out(self.fitter)
      npreds, nstr = errors.shape
      prediction = errors[ numpy.arange(errors.shape[0]), numpy.arange(errors.shape[0]) ]
      return numpy.average(prediction*prediction)
    else:
      errors = ce.leave_many_out(self.fitter, self.sets)
      npreds, nstr = errors.shape
      prediction = errors[ [[(j,i) for i in self.sets[j]] for j in len(self.sets)] ]
      return numpy.average(prediction*prediction)



  def __call__(self, indiv):
    from math import fabs

    result = 0e0
    for r in indiv.genes - self.target:
      result += fabs(r)
    return float(result) / float(len(indiv.genes))

    

