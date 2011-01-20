""" Bandstructure plotting tools """
__docformat__  = 'restructuredtext en'
__all__ = ['band_structure', 'Extract']

from ..opt.decorators import broadcast_result, make_cached
from ._extract import MassExtract

def band_structure(escan, structure, kpoints, density = None, outdir=None, comm=None,\
                   do_relax=None, pools = 1, nbkpoints = None, **kwargs):
  """ Returns 
  
      :Parameters:
        escan : `lada.escan.Escan`
          Functional wrapping nanopse's ESCAN. 
        structure : `lada.crystal.Structure`
          object describing the structure for which to compute a
          band-structure.
        kpoints : sequence of 2-tuples
          Each two tuple is a starting k-point and an end k-point. The k-points
          should be given in cartesian units within the ideal, undistorted
          lattice. They will deformed to the fit into the relaxed structure.
          However, translational invariance is not applied (e.g. kpoints are
          not refolded).
        density : float 
          Number of kpoints per reciprocal length unit.
        kwargs 
          Any parameters appropriate to `lada.escan.Escan`.

      :return: sequence of (x, kpoint, eigenvalues), where
        - x is the abscissa for the brand-structure plot.
        - kpoint is actual deformed kpoint at which computation was performed.
        - eigenvalues is a numpy array of eigenvalues.
  """
  from os import getcwd
  from os.path import join, expanduser, abspath, exists
  from shutil import copyfile
  from boost.mpi import world, all_gather, broadcast
  from numpy.linalg import norm
  from numpy import abs, sum
  from ..crystal import deform_kpoint
  from ..opt import RelativeDirectory

  # check/correct input arguments
  assert nbkpoints == None or density == None, ValueError("Choose either density or nbkpoints")
  assert nbkpoints != None or density != None, ValueError("Choose either density or nbkpoints")
  assert "do_genpot" not in kwargs,\
         ValueError("\"do_genpot\" is not an admissible argument of band_structure.")
  assert "do_escan" not in kwargs,\
         ValueError("\"do_escan\" is not an admissible argument of band_structure.")
  outdir = RelativeDirectory(outdir if outdir != None else getcwd()).path
  outdir_calc = join(outdir, 'calculations')
  if comm == None: comm = world
  if pools > comm.size: pools = comm.size
  vffrun = kwargs.pop("vffrun", escan.vffrun)
  genpotrun = kwargs.pop("genpotrun", escan.genpotrun)

  # first computes vff and genpot unless given.
  if genpotrun == None or vffrun == None: 
    vffout = escan( structure, outdir=outdir_calc, do_escan=False, genpotrun=genpotrun,\
                    vffrun=vffrun, comm = comm, **kwargs )
    if genpotrun == None: genpotrun = vffout
    if vffrun == None: vffrun = vffout
  
  # two functions required to continue.
  input, relaxed = structure.cell.copy(), vffout.structure.cell.copy()
  dont_deform_kpoint = vffout.functional._dont_deform_kpoint
  def _get_kpoint(_kpoint):
    """ Deforms kpoint to new lattice, if required. """
    if dont_deform_kpoint: return _kpoint
    if sum(abs(input-relaxed)) < 1e-11: return _kpoint
    return deform_kpoint(_kpoint, input, relaxed)
  def _line(start, end):
    """ Generator for creating points between two kpoints. """
    from numpy.linalg import norm

    distance = norm(_get_kpoint(end - start))
    if nbkpoints != None: nbkpt = nbkpoints
    else: nbkpt = int(max(1, float(density) * distance - 1))
    stepsize = 1e0/float(nbkpt)
    _kpoints = [ float(i) * stepsize for i in range(1, nbkpt+1) ]
    for k in _kpoints: yield start + k * (end-start) 

  def _lines(endpoints):
    """ Generator for creating segments. """
    from numpy.linalg import norm
    assert len(endpoints) > 0, ValueError
    assert len(endpoints[0]) == 2, ValueError
    pos = 0
    yield pos, endpoints[0][0]
    for start, end in endpoints:
      last = start.copy()
      for _kpoint in _line(start, end):
        pos += norm(_get_kpoint(_kpoint-last))
        last = _kpoint.copy()
        yield pos, _kpoint


  # splits local communicator.
  color = comm.rank % pools
  local_comm = comm.split(color)
  # then computes different kpoints.
  results = []
  for i, (x, kpoint) in enumerate(_lines(kpoints)):
    # separates jobs into pools.
    if i % pools != color: continue

    # sets directory.
    directory = join(outdir_calc, "%i-%s" % (i, kpoint))
    # actually computes stuff.
    out = escan( structure, outdir=directory, kpoint=kpoint, vffrun=vffrun,\
                 genpotrun=genpotrun, do_escan=True, comm = local_comm, **kwargs )
    # saves stuff
    eigenvalues = out.eigenvalues.copy()
    eigenvalues.sort()
    results.append( (x, kpoint, eigenvalues) )

  if comm.size > 1 and pools > 1: # gathers and orders results.
    head_comm = comm.split(0 if local_comm.rank == 0 else 1)
    if local_comm.rank == 0:
      results = all_gather(head_comm, results)
      def comp(a,b):
        if a[0] < b[0]: return -1
        if a[0] ==  b[0]: return 0
        return 1
      results = sorted((j for i in results for j in i), comp)
      broadcast(local_comm, results, 0)
    else: results = broadcast(local_comm, None, 0) 
  return results

class Extract(MassExtract):
  """ Extraction class for band-structures. """

  def __iter_alljobs__(self):
    """ Goes through all calculations and orders them. """
    from re import compile
    from operator import itemgetter

    alls = []
    regex = compile("/calculations/(\d+)-\[\s*\S+\s+\S+\s+\S+\s*\]")
    for key, result in MassExtract.__iter_alljobs__(self):
      found = regex.match(key)
      if found != None: alls.append((key, result, i))
    for key, value, i in result: yield key, value

  @property
  @make_cached
  def success(self): 
    """ Checks for success of jobs. """
    from os.path import exists
    if not exists(self.rootdir): return False
    if len(self.items()) == 0: return False
    for name, job in self.iteritems(): 
      if not job.success: return False
    return True

  @property
  @make_cached
  def vff(self):
    """ Vff extraction object. """
    from os.path import join
    return self.Extract(join(self.rootdir, 'calculations'))

  @property
  @make_cached
  def kpoints(self):
    """ kpoints used to compute band-structure. """
    from numpy import zeros, array
    from re import compile
    regex = compile("/calculations/(\d+)-\[\s*(\S+)\s+(\S+)\s+(\S+)\s*\]")
    result = zeros((len(self.jobs), 3), dtype='float64')
    for i, name in enumerate(self.iterkeys()):
      found = regex.search(name)
      result[i,:] = array([found.group(2), found.group(3), found.group(4)])
    return result

  @property
  def vbm(self): 
    """ Returns energy at vbm. """
    from numpy import array, max
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.vff.structure)
    units = self.eigenvalues.itervalues().next().units
    return max(array(self.eigenvalues.values())[:, nbe-2:nbe]) * units

  @property
  def cbm(self): 
    """ Returns energy at vbm. """
    from numpy import array, min
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.vff.structure)
    units = self.eigenvalues.itervalues().next().units
    return min(array(self.eigenvalues.values())[:, nbe:nbe+2]) * units

  @property 
  def directness(self):
    """ Difference in energy between the CBM at Gamma and the LUMO. """
    from numpy.linalg import norm
    from ..crystal import nb_valence_states
    lumo = self.cbm
    gamma = min((job for job in self.values()), key=lambda x: norm(x.functional.kpoint))
    if norm(gamma.functional.kpoint) > 1e-6: raise RuntimeError("Gamma point not found.")
    nbe = nb_valence_states(self.vff.structure)
    cbm = min(gamma.eigenvalues[nbe], gamma.eigenvalues[nbe+1])
    return cbm - lumo
    
try: import matplotlib.pyplot as plt 
except: 
  def plot_bands(extractor, **kwargs):
    """ Plots band-structure. """
    raise ImportError("Cannot use plot_bands without matplotlib. """)
else:
  def plot_bands(extractor, tolerance=1e-6, **kwargs):
    """ Tries and plots band-structure. """
    from numpy import dot, array, min, max
    from numpy.linalg import norm

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
    y = array(extractor.eigenvalues.values())

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

  Extract.plot_bands = plot_bands



