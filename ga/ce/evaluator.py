""" Holds evaluation functions for cluster expansion. """
__docformat__ = "restructuredtext en"


class Evaluator(object):
  """ Evaluator class for cluster expansion. """
  def __init__(self, lattice, structures, energies=None, lmo=0.33333, nsets=5, **keywords):
    """ Creates and evaluator for cluster expansion.  
        
        :Parameters:
          lattice 
            Back-bone lattice for the generalized Ising model.
          structures
            List of structures forming the fitting and training sets.
          energies 
            List of energies for the structures. Alternatively, the energies
            can be taken from the corresponding attibute of each structure in
            the list.
          lmo : float
            ratio of training vs validation structures for each set.
          nsets : int
            number of training vs validation sets.
    """
    from re import compile, match
    from random import shuffle
    from numpy import array
    from ...ce import cluster_factory

    super(Evaluator, self).__init__()

    self.lattice = lattice
    """ Back-bone lattice of the Ising model. """
    
    # creates clusters.
    args = {}
    key_regex = compile("B(\d+)")
    for key, value in keywords.iteritems(): 
      if key == "J0": args[key] = value
      elif key == "J1": args[key] = value
      elif match(key_regex, key) != None: args[key] = value
    for key in args.iterkeys(): del keywords[key]
    self._clusters = cluster_factory(self.lattice, **args) 
    """ List of clusters. """
    self._args = args
    """ Holds arguments to create clusters. """

    self._structures = structures
    """ List of structures. """
    self._pis = array([self.clusters.pis(s, self.lattice) for s in self.structures], dtype="float64")
    """ List of pis for input structures. """

    # creates or assigns energies as needed.
    self._energies = energies
    """ List of energies for each structure. """
    if self._energies == None: self._energies = array([s.energy for e in self.structures])
    assert len(self.energies) == len(self.structures),\
           ValueError("Inequivalent number of structures and energies on input.")

    # creates sets.
    assert lmo > 0e0 and lmo < 1e0, ValueError("Leave-many-out ratio should be between 0 and 1.")
    crosets = []
    n = min(max(1, int(lmo * float(len(self.structures)) + 0.5)), len(self.structures))
    indices = []
    for i in range(nsets):
      if len(indices) == 0:
        indices = range(len(self.structures))
        shuffle(indices)
      if len(indices) >= n:
        crosets.append(indices[:n])
        indices = indices[n:]
      if len(indices) < n: 
        crosets.append(indices)
        u = n - len(crosets[-1])
        indices = range(len(self.structures))
        shuffle(indices)
        crosets[-1].extend(indices[:u])
        indices = indices[u:]
    self.crosets = array(crosets)
    """ List of cross-validation sets. """

      
  @property
  def clusters(self):
    """ List of MBCE clusters. """
    return self._clusters

  @property
  def structures(self):
    """ List of structures over which to perform fit.. """
    return self._structures
    
  @property
  def pis(self):
    """ List of pis for known list of structures. """
    return self._pis

  @property
  def energies(self):
    """ List of energies for known list of structures. """
    return self._energies
  @energies.setter
  def energies(self, values):
    """ List of energies for known list of structures. """
    from numpy import array
    values = array(values).copy()
    assert len(values) == len(self.structures),\
           ValueError("There should be as many energies as structure.")
    self._energies = values

  @property
  def crosets(self):
    """ List of sets for which to perform cross-validation. """
    return self._crosets
  @crosets.setter
  def crosets(self, value):
    """ Sets list of cross-validation sets. """
    from numpy import argmin, argmax, array
    assert argmin(value) >= 0, ValueError("Cross-validation sets cannot contain negative indices.")
    assert argmax(value) < len(self.structures),\
           ValueError("Cross-validation sets cannot contain indices "\
                      "larger than the number of structures.")
    self._crosets = array(value).copy()
    self._fitsets = array([ [i for i in range(len(self.structures)) if i not in croset] \
                            for croset in self.crosets ])
  
  @property
  def fitsets(self):
    """ List of sets for which to perform fits. """
    return self._fitsets


  def run(self, genes):
    """ Computes cross validations for given set of clusters. """
    from numpy import sum, dot, array
    from numpy.linalg import lstsq

    # creates boolean array from genes.
    genes = array([i in genes for i in xrange(self.pis.shape[1])])
    # selects pis for given set of clusters.
    pis = self.pis[:,genes]

    # loops over cross-validation sets.
    scores = []
    for croset, fitset in zip(self.crosets, self.fitsets): 
      # pis, energies in cross-validation set.
      cropis = pis[croset, :]
      croene = self.energies[croset]
      # pis, energies in fitting set.
      fitpis = pis[fitset,:]
      fitene = self.energies[fitset]
      # performs fitting.
      interactions = lstsq(fitpis, fitene)[0]
      scores.append(sum((dot(cropis, interactions) - croene)**2))
      
    return scores

  def __call__(self, indiv, comm=None, outdir=None, **kwargs):
    """ Evaluates cross validation scores for current gene-set. """
    from numpy import mean
    this = self.copy(**kwargs)
    indiv.cvscores = this.run(indiv.genes)
    indiv.fitness = mean(indiv.cvscores)
    return indiv.fitness

  def __repr__(self):
    """ Representation of this structure. """
    string += repr(self.lattice)
    string += "structures = []\n"
    for structure in self.structure:
      string += repr(structure)
      string += "structures.append(structure)\n"
    string += "evaluator = {0.__class__.__name__}(lattice, structures, ".format(self)
    for key, value in self._args: string += "{0}={1}".format(key, value)
    string += ")\n"
    string += "evaluator.crosets = {0}\n".format(repr(list(self.crosets)))
    string += "evaluator.energies = {0}\n".format(repr(list(self.energies)))
    return string

  def copy(self, **kwargs):
    """ Performs deep-copy of object. 

				:Param kwargs: keyword arguments will overide the corresponding
                 attribute. The value of the keyword argument is *not* deepcopied. 
                 The attribute must exist. Attributes cannot be created here. 
    """
    from copy import deepcopy
    result = deepcopy(self)
    noadd = kwargs.pop('noadd', None)
    for key, value in kwargs.iteritems():
      if not hasattr(result, key):
        if noadd == None:
          raise ValueError( "Attribute {0} does not exist "\
                            "and will not be created in copy.".format(key))
        elif noadd == True: continue
      setattr(result, key, value)
    return result

class LocalSearchEvaluator(Evaluator):
  """ Performs simple minded local search. """
  def __init__( self, lattice, structures, taboos=None, energies=None, lmo=0.33333, nsets=5,
                maxiter=-1, maxeval=-1, exclude=None, **keywords):
    """ Initializes the local search functional.
    
        :Parameters:
          lattice 
            Back-bone lattice for the generalized Ising model.
          structures
            List of structures forming the fitting and training sets.
          energies 
            List of energies for the structures. Alternatively, the energies
            can be taken from the corresponding attibute of each structure in
            the list.
          lmo : float
            ratio of training vs validation structures for each set.
          nsets : int
            number of training vs validation sets.
          maxiter : int
            Maximum number of iterations.
          maxeval : int
            Maximum number of evaluations.
          exclude : sequence or None
            Clusters to exclude from optimization.
    """
    super(LocalSearchEvaluator, self).__init__(lattice, structures, energies, lmo, nsets, **keywords)
    self.maxiter = maxiter
    """ Maximum number of iterations. """
    self.maxeval = maxeval
    """ Maximum number of evaluations. """
    self.exclude = set(exclude) if exclude != None else set()
    """ Clusters/bits to exclude from optimization. """

  def __call__(self, indiv, comm=None, outdir=None, verbose=False, **kwargs):
    """ Performs a simple minded local search. """
    # case without optimization.

    if self.maxiter == 0:
      if verbose: print "no local search."
      return super(LocalSearchEvaluator, self).__call__(indiv, comm=comm, outdir=outdir, **kwargs)

    from random import shuffle
    from copy import deepcopy
    indiv = deepcopy(indiv)

    if not hasattr(indiv, "fitness"):
      super(LocalSearchEvaluator, self).__call__(indiv, comm=comm, outdir=outdir, **kwargs)
    iteration, evaluation = 0, 0
    indiv.genes = set(indiv.genes)
    if verbose:
      print "Local search:\n  maxiter = {0}\n  maxeval {1}\n start {2} -- {3}\n"\
            .format(self.maxiter, self.maxeval, indiv.genes, indiv.fitness)

    indices = [u for u in xrange(indiv.maxsize) if u not in self.exclude]
    while     (self.maxiter < 0 or iteration < self.maxiter) \
          and (self.maxeval < 0 or evaluation < self.maxeval):

      if verbose: print "** Iteration {0} **".format(iteration)

      changed = False
      shuffle(indices)
      current = indiv.fitness
      for i in indices:
        other = deepcopy(indiv)
        if i in other.genes: other.genes.remove(i)
        else: other.genes.add(i)
        if self.darwin.taboo(other): continue
        super(LocalSearchEvaluator, self).__call__(other, comm=comm, outdir=outdir, **kwargs)
        self.darwin.taboo(other) # update with fitness for later analysis.
        evaluation += 1
        if self.darwin.comparison(indiv, other) == 1: 
          indiv = other
          changed = True
        if self.maxeval > 0 and evaluation < self.maxeval: break
      if verbose: print "change {0}, current {1}.".format(indiv.fitness - current, indiv.fitness)

      if not changed: break
      iteration += 1

    indiv.genes = sorted(indiv.genes)
    if verbose: print "end local search: {0} {1}.\n".format(indiv.genes, indiv.fitness)
    return indiv.fitness

  def __repr__(self):
    """ Throws! Cannot represent itself. """
    raise NotImplementedError("Evaluation object cannot represent itself. """)
