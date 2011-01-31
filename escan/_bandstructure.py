""" Bandstructure plotting tools """
__docformat__  = 'restructuredtext en'
__all__ = ['band_structure', 'Extract']

from ..opt.decorators import broadcast_result, make_cached
from ._extract import MassExtract
from .kpoints import KPoints, _reduced_grids_factory

class BPoints(KPoints):
  """ *Grid* of kpoints for bands-structure calculation. """
  def __init__(self, lines=None, density=10, relax=True, mink=3):
    KPoints.__init__(self)
    self.lines = lines
    """ Lines of k-points for which to perform calculations. """
    if self.lines == None: self.lines = []
    self.density = density
    """ Requested density of kpoints. """
    self.mink = mink
    """ Minimum number of kpoints per line. """
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """

  def _mnk(self, input, output):
    """ Yields lines of k-points with appropriate density. """
    from numpy import any, abs, dot
    from numpy.linalg import norm, inv

    if len(self.lines) == 0: return
    if self.relax: deformation = dot(inv(output.cell.T), input.cell.T)
    
    last_end = self.lines[0][0] + 1e0
    for start, end in self.lines:
      if any(abs(last_end-start) > 1e-12):
        yield 1e0, (dot(deformation, start) if self.relax else start)
      last_end = end
      for kpoint in self._line(start, end):
        yield 1e0, (dot(deformation, kpoint) if self.relax else kpoint)

  def _line(self, start, end):
    """ Yields kpoints on a line. """
    from numpy.linalg import norm
    distance = norm(end-start)
    nbkpt = max(self.mink, int(float(self.density) * distance - 1))
    stepsize = 1e0 / float(nbkpt)
    for i in range(1, nbkpt+1):
      yield start + float(i) * stepsize * (end-start)

  def __iadd__(self, values):
    """ Adds segment to band-structure line. """
    from numpy import array
    assert hasattr(values, '__iter__'), ValueError('Value should be a sequence.')
    values = [v for v in values]
    if len(values) == 2 and len(values[0]) == 3 and len(values[1]) == 3: values = [values]
    for start, end in values:
      assert hasattr(start, '__iter__'), ValueError('Values should be two tuples of 3d-vectors.')
      start = array(start, dtype='float64')
      assert hasattr(end, '__iter__'), ValueError('Values should be two tuples of 3d-vectors.')
      end = array(end, dtype='float64')
      self.lines.append((start, end))
    return self

  def __add__(self, other):
    """ Returns new instance with summed segment. """
    result = self.__class__(self.lines, self.density, self.relax, self.mink)
    result += other.lines if hasattr(other, 'lines') else other
    return result

  def __repr__(self):
    """ Returns representation of this object. """
    compare = BPoints()
    string = '{0.__class__.__name__}(None, {0.density}, {0.relax}, {0.mink})'.format(self)
    for start, end in self.lines:
      string += '+([{0[0]},{0[1]},{0[2]}], [{1[0]},{1[1]},{1[2]}])'.format(start, end)
    return string
    

ReducedBPoints = _reduced_grids_factory('ReducedBPoints', BPoints)

try: import matplotlib.pyplot as plt 
except: 
  def plot_bands(extractor, **kwargs):
    """ Plots band-structure. """
    raise ImportError("Cannot use plot_bands without matplotlib. """)
  def plot_alloybands(extractor, multicell, tolerance=1e-6, **kwargs):
    """ Plots alloy band-structure using the majority representation. """
    raise ImportError("Cannot use plot_alloybands without matplotlib. """)
else:
  def plot_bands(extractor, **kwargs):
    """ Tries and plots band-structure. """
    from numpy import dot, array, min, max
    from numpy.linalg import norm

    old = extractor.unreduce
    extractor.unreduce = True

    bandcolor = kwargs.pop('bandcolor', 'blue')
    edgecolor = kwargs.pop('edgecolor', 'red')
    edgestyle = kwargs.pop('edgestyle', '--')

    # first finds breaking point.
    kpoints = extractor.kpoints
    delta = kpoints[1:] - kpoints[:-1]
    norms = [norm(delta[i,:]) for i in range(delta.shape[0])]
    bk = []
    for i, d in enumerate(norms[1:]):
      if abs(norms[i]-d) > 1e-6: bk.append(i+1)

    # then plot bands.
    x = array([sum(norms[:i]) for i in range(len(norms)+1)])
    y = array(extractor.eigenvalues)

    # then line markers.
    plt.plot(x, y, color=bandcolor, **kwargs)
    for i in bk: plt.axvline(x[i], color='black', **kwargs)

    # then plot vbm and cbm.
    kwargs.pop('linestyle', None) 
    plt.axhline(extractor.vbm, color=edgecolor, linestyle=edgestyle, **kwargs)
    plt.axhline(extractor.cbm, color=edgecolor, linestyle=edgestyle, **kwargs)



    plt.xlim((x[0], x[-1]))
    ylims = min(y) - (max(y) - min(y))*0.05, max(y) + (max(y) - min(y))*0.05
    plt.ylim(ylims)

    extractor.unreduce = old

  def plot_alloybands(extractor, multicell, tolerance=1e-6, **kwargs):
    """ Plots alloy band-structure using the majority representation. """
    from numpy import dot, array, min, max
    from numpy.linalg import norm
    from . import majority_representation
    from pickle import load, dump

    edgecolor = kwargs.pop('edgecolor', 'red')
    edgestyle = kwargs.pop('edgestyle', '--')

    # first finds breaking point.
    istr = extractor.input_structure
    ostr = extractor.structure
    kpoints = array([u[1] for u in extractor.functional.kpoints.unreduced(istr, ostr)])
    delta = kpoints[1:] - kpoints[:-1]
    norms = [norm(delta[i,:]) for i in range(delta.shape[0])]
    bk = []
    for i, d in enumerate(norms[1:]):
      if abs(norms[i]-d) > 1e-6: bk.append(i+1)

    # then plot bands.
    xvalues = array([sum(norms[:i]) for i in range(len(norms)+1)])
    args = [[],[],[]]
    maj = majority_representation(extractor, multicell)
    for values, x in zip(maj, xvalues):
      for v in values: 
        args[0].append(x)
        args[1].append(v[0])
        args[2].append(v[1])

    # then line markers.
    plt.hexbin(*args, **kwargs)
    x = array([sum(norms[:i]) for i in range(len(norms)+1)])
    for i in bk: plt.axvline(x[i], color='black', **kwargs)

    kwargs.pop('linestyle', None) 
    plt.axhline(extractor.vbm, color=edgecolor, linestyle=edgestyle, **kwargs)
    plt.axhline(extractor.cbm, color=edgecolor, linestyle=edgestyle, **kwargs)

    y = array(extractor.eigenvalues)
    plt.xlim((x[0], x[-1]))
    ylims = min(y) - (max(y) - min(y))*0.05, max(y) + (max(y) - min(y))*0.05
    plt.ylim(ylims)
