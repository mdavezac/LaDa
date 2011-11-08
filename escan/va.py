""" Allows for construction of Virtual Atom Derivatives. """
from boost import mpi

class VA(object):

  vff_inputfile = "atom_input." + str( mpi.world.rank )

  def __init__(self, functional, vff, outsize, dir = None, order = 2):
    from copy import deepcopy

    self.functional = deepcopy(functional)
    self.vff = deepcopy(vff)
    self.outsize = outsize
    self._size = None
    self._which = None
    self.dir = dir
    self.order = 2


  def init(self, structure, dir = None, order = 2):
    """ Initializes for a given structure. """
    from copy import deepcopy
    import os
    import tempfile

    if dir is not None: self.dir = dir
    if order != self.order: self.order = order

    self.vff.structure = structure
    lattice = structure.lattice
    self._size = 0
    self._which = []
    for i, atom in enumerate(structure.atoms):
      assert atom.site >= 0 and atom.site < len(lattice.sites)
      if len(lattice.sites[atom.site]) > 1:
        self._size += 1
        self._which.append(i)

    if self.dir is None:
      tempfile.tempdir = self.directory
      self.dir = tempfile.mkdtemp()

    if not os.path.exists(self.dir): 
      os.path.makedirs(self.dir)

  def __len__(self): return self._size

  def _setx(self, x):
    from math import fabs

    result = True
    for i in which:
      self.vff.structure.atoms[i].type = x
      if fabs(x) > 1e-8: result = False
    return result

  def _print_escan_input():
    from os.path import join as make_path
    from lada.physics import a0

    file = open(make_path(self.dir, self.vff_inputfile), "w")
    for i in xrange(3):
      for j in xrange(3):
        print >>file, structure.cell(i,j) * structure.scale / a0("A"),
      for j in xrange(3):
        print >>file, structure.cell(i,j)
      print >>file



  def __call__(self, x):
    """ Returns value of function at x. """
    import numpy as np
    from os.path import join as make_path

    assert self._size is not None

    run_vff = self._setx(x)

    if run_vff:
      vff.evaluate()
      vff.print_escan_input( make_path(self.dir, self.vff_inputfile) )
    else: self._print_escan_input()

    return self.functional()


def add_va_gradient(escan, vff, name="gradient", outsize=1 ):
  """ Adds an VA gradient bound-method to self. """
  import type

  def gradient(self, structure, dir=None, order = 2):
    """ A virtual atom gradient method for escan and friends.
        structure is the structure for which to compute VA gradients.
    """
    import numpy as np
    import os
    import tempfile

    # gets size of return from structure.
    result = np.matrix( np.zeros((size, outsize), dtype="float64") )

    # computes order 0 (unless path given).
    self.directory = dir + "_a0"
    a0 = self( structure )

    for i in xrange(size):
      values = np.zeros( (outsize, order+1), dtype="float64" )









  setattr(self, name, types.MethodType(gradient, self, self.__class__))
