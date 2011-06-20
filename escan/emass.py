""" Functional to compute effective mass. """
__docformat__ = "restructuredtext en"
__all__ = ["Functional", "Extract"]
from ..opt.decorators import make_cached
from ._bandgap import Functional as Bandgap
from .kescan import Extract as KExtract

class Extract(KExtract):
  """ Extraction object for effective mass functional. """
  def __init__(self, directory=None, comm=None, bandgap=None, **kwargs):
    """ Initializes the extraction object. """
    super(Extract, self).__init__(directory, comm, **kwargs)
    if bandgap != None: self._extract_bg = bandgap
  
  @property
  def type(self): 
    """ Whether electronic or hole effective masses are computed. """
    return self.functional.type

  @property
  def direction(self): 
    """ Direction for which effective mass is computed. """
    return self.functional.direction

  @property
  def center(self): 
    """ K-point at which effective mass is computed. """
    return self.functional.center 

  @property
  def extract_bg(self, *args, **kwargs):
    """ Computes dipole element between vbm and cbm. """
    if hasattr(self, "_extract_bg"): return self._extract_bg
    from . import extract_bg as ebg
    e = ebg(self.directory, self.comm)
    assert e.success, RuntimeError("Could not extract bandgap from calculations.")
    return e

  @property 
  def cbm(self):
    """ Conduction Band Minimum at kpoint for which effective mass is computed. """
    return self.extract_bg.cbm 

  @property 
  def vbm(self):
    """ Valence Band Minimum at kpoint for which effective mass is computed. """
    return self.extract_bg.vbm 

  @property 
  def bandgap(self):
    """ Band-gap at kpoint for which effective mass is computed. """
    return self.extract_bg.bandgap
  
  def dipole(self, *args, **kwargs):
    """ Computes dipole element between vbm and cbm. """
    return self.extract_bg.dipole(*args, **kwargs)

  def oscillator_strength(self, *args, **kwargs):
    """ Computes dipole element between vbm and cbm. """
    return self.extract_bg.dipole(*args, **kwargs)

  @property
  @make_cached
  def derivatives(self): 
    """ Effective mass in given direction. 
    
        Returned in units of the electron's mass at rest. Sign difference
        between electron and holes is accounted for.
    """
    from numpy import sort, array
    from numpy.linalg import lstsq as np_lstsq
    from quantities import hartree
    lstsq = self.functional.lstsq
    if lstsq == None: lstsq = np_lstsq

    # constructs matrix of eigenvalues with cbm or vbm only.
    eigs = []
    eigenvalues = KExtract.eigenvalues.__get__(self).rescale(hartree)
    if self.type == "e": 
      for e in eigenvalues:
        eigs.append( e[e > (self.vbm + 0.5 * self.bandgap)] )
    elif self.type == "h":
      for e in eigenvalues:
        eigs.append( e[e < (self.vbm + 0.5 * self.bandgap)] )
    else: raise RuntimeError("type ({0}) is neither \"e\" nor \"h\".")
    m = min([len(e) for e in eigs])
    eigs = array([ e[:m] for e in eigs ]) * hartree
    measurements = sort(eigs, axis=1) 
    for k in self.functional.kpoints(self.input_structure, self.structure): continue
    return lstsq(self.functional.kpoints.parameters, measurements)[0]

  @property
  def eigenvalues(self): 
    """ Eigenvalues at k-point where effective masses are computed. 
    
        These eigenvalues are the result of the least-square-fit.
    """
    from quantities import eV, hartree
    result = (self.derivatives[0,:] * hartree).rescale(eV)
    return result if self.type == "e" else result[::-1]

  @property
  def mass(self): 
    """ Effective mass in given direction. 
    
        Returned in units of the electron's mass at rest. Sign difference
        between electron and holes is accounted for.
    """
    assert self.type == "e" or self.type == "h",\
           RuntimeError("Unknown type {0}.".format(self.type))
    return 1e0/self.derivatives[2,:] if self.type == "e" else -1e0/self.derivatives[2,::-1]

class Functional(Bandgap):
  """ Effective mass functional. 
  
      The effective mass is computed for a given set of bands, for the VBM, or
      the CBM, in a  specific direction. 
  """
  Extract = Extract
  """ Extraction object for the effective-mass functional. """
  def __init__( self, direction=(0,0,1), nbpoints=None, stepsize=1e-2, \
                center=None, lstsq=None, **kwargs ):
    """ Computes effective mass for a given direction.
    
        :Parameters:
          type : "e" or "h"
            Whether to compute electronic or effective mass.
          direction : 3-tuple 
            direction for which to compute effective mass.
          nbpoints : int
            Number of points with wich to compute derivatives.
            Should be at least order + 1. Default = order + 1. 
          stepsize : float
            Distance between interpolation points. Default = 1e-2.
            Units of ``2|pi|/a``, with ``a=structure.scale``.
          center : 3-tuple
            k-point where to take derivative. Units of ``2|pi|/a``, with
            ``a=structure.scale``.
          lstsq 
            Linear least square method. The first two parameters should
            be same as numpy.linalg.lstsq. Other parameters can be passed as extra
            parameters. Defaults to numpy.linalg.lstsq.
          kwargs 
            Extra parameters which are passed on first to escan (if escan
            object has an attribute with the same name), then to the linear least
            square fit method. Note that this will *not* change the external escan
            object.  This function is stateless. 
    
        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    self.type      = "e"
    """ Whethe to compute electronic or hole effective masses. """
    self.direction = direction
    """ Direction for which to compute effective mass. """
    self.nbpoints = nbpoints
    """ Number of points with which to perform least-square fit. Defaults to 3. """
    self.stepsize = stepsize
    """ Distance between interpolation points. Default = 1e-2.

        Units of ``2|pi|/a``, with ``a=structure.scale``.

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    self.center = center
    """ k-point for which to compute effective mass. """
    self.lstsq = lstsq
    """ Least square fit method to use when computing effective mass. 
    
        If None, will default ot numpy.linalg.lstsq. Otherwise, it should be a
        pickleable-callable, or suffer the consequences.
    """
    super(Functional, self).__init__(**kwargs)

  def __call__(self, structure, outdir=None, comm=None, bandgap=None, **kwargs):
    """ Computes effective mass.

        :Parameters: 
          structure : `crystal.Structure`
            Structure for which to compute effective mass.
          outdir : str or None
            Directory where to save results. If None, will create results are
            stored in current directory.
          comm : None or `mpi.Communicator`
            MPI Communicator holding processes with which to perform
            calculations.
          bandgap 
            If given, should be the return from a bandgap calculation.
        Parameters are passed on to `dervatives.reciprocal` method.
    """
    from os import symlink
    from os.path import join, exists, relpath, dirname
    from ..opt import copyfile
    from .derivatives import reciprocal
    if '_computing' in self.__dict__:
      return super(Functional, self).__call__(structure, outdir, comm, **kwargs)

    type      = kwargs.pop('type', self.type)
    direction = kwargs.pop('direction', self.direction)
    nbpoints  = kwargs.pop('nbpoints', self.nbpoints)
    stepsize  = kwargs.pop('stepsize', self.stepsize)
    center    = kwargs.pop('center', self.center)
    lstsq     = kwargs.pop('lstsq', self.lstsq)

    # copy functional with current type
    this = self.copy(type=type)

    # first computes bandgap.
    if bandgap != None:  # in this case, just links to previous calculations.
      assert bandgap.success,\
             RuntimeError( "Input bandgap calculations at {0} were not successfull."\
                           .format(bandgap.directory) )
      result = self.Extract(outdir, comm, unreduce=True, bandgap=bandgap)
      directory = bandgap.directory
      if bandgap.is_ae: directory = dirname(directory)
      for file in bandgap.iterfiles(): 
        name = join(result.directory, relpath(file, bandgap.directory))
        print file, name, relpath(file, bandgap.directory), bandgap.__class__.__module__
        if not exists(name): # tries and links to earliest directory.
          alldirs = relpath(file, bandgap.directory).split('/')
          src, dest, done = bandgap.directory, result.directory, False
          for now in alldirs:
            src, dest = join(src, now), join(dest, now)
            if not exists(dest):
              copyfile(src, dest, symlink=True)
              done = True
              break
    assert False

    bandgap = super(Functional, this).__call__(structure, outdir, comm, **kwargs)

    assert bandgap.success, RuntimeError("Could not compute bandgap.")
    
    # then figures out appropriate reference.
    if type == "e":   eref = bandgap.cbm - bandgap.bandgap / 4e0
    elif type == "h": eref = bandgap.vbm + bandgap.bandgap / 4e0
    
    # performs calculation.
    reciprocal(this, structure, outdir, direction=direction,
               nbpoints=nbpoints, stepsize=stepsize, center=center, 
               lstsq=lstsq, eref=eref, order=2, **kwargs)

    # creates extraction object.
    result = self.Extract(outdir, comm, unreduce=True, bandgap=bandgap)

    # makes sure bandgap directory exists.
    if bandgap != None: 
      if bandgap.is_ae:
        if not exists(join(result.directory, "AE")):
          symlink(join(bandgap.directory, "AE"), join(result.directory, "AE"))
      else: 
        if not exists(join(result.directory, "VBM")):
          symlink(join(bandgap.directory, "VBM"), join(result.directory, "VBM"))
        if not exists(join(result.directory, "CBM")):
          symlink(join(bandgap.directory, "CBM"), join(result.directory, "CBM"))
         
        if not exists(join(result.directory, "vff_out")):
          symlink(join(bandgap.directory, "vff_out"), join(result.directory, "CBM"))

    # Effective mass extractor.
    return self.Extract(outdir, comm, unreduce=True, bandgap=bandgap)

