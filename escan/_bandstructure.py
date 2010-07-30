#! python
""" Bandstructure plotting tools """
__docformat__  = 'restructuredtext en'

def band_structure(escan, structure, kpoints, density, outdir=None, comm=None,\
                   do_relax=None, pools = 1, **kwargs):
  """ Returns eigenvalues for plotting bandstructure. 
  
      :Parameters:
      - `escan`: an `lada.escan.Escan` functional wrapping nanopse's pescan. 
      - `structure`: an `lada.crystal.Structure` object describing the
        structure for which to compute a band-structure.
      - `kpoints`: a sequence of 2-tuples. Each two tuple is a starting k-point and
        an end k-point. The k-points should be given in cartesian units within
        the ideal, undistorted lattice. They will deformed to the fit into the
        relaxed structure. However, translational invariance is not applied
        (e.g. kpoints are not refolded).
      - `density`: a float giving the number of kpoints per reciprocal length
        unit.
      - **kwargs: Any parameters appropriate to `lada.escan.Escan`.
  """
  from os import getcwd
  from os.path import join, expanduser, abspath, exists
  from shutil import copyfile
  from boost.mpi import world, all_gather, broadcast
  from numpy.linalg import norm
  from numpy import abs, sum
  from ..crystal import deform_kpoint

  # check/correct input arguments
  assert "do_genpot" not in kwargs,\
         ValueError("\"do_genpot\" is not an admissible argument of band_structure.")
  assert "do_escan" not in kwargs,\
         ValueError("\"do_escan\" is not an admissible argument of band_structure.")
  outdir = abspath(expanduser(outdir)) if outdir != None else getcwd()
  outdir = join(outdir, "band_structure")
  if comm == None: comm = world
  if pools > comm.size: pools = comm.size
  vffrun = kwargs.pop("vffrun", escan.vffrun)
  genpotrun = kwargs.pop("genpotrun", escan.genpotrun)

  # first computes vff and genpot unless given.
  if genpotrun == None or vffrun == None: 
    vffout = escan( structure, outdir=outdir, do_escan=False, genpotrun=genpotrun,\
                    vffrun=vffrun, comm = comm, **kwargs )
    if genpotrun == None: genpotrun = vffout
    if vffrun == None: vffrun = vffout
  
  # two functions required to continue.
  input, relaxed = structure.cell.copy(), vffout.structure.cell.copy()
  dont_deform_kpoint = vffout.escan._dont_deform_kpoint
  def _get_kpoint(_kpoint):
    """ Deforms kpoint to new lattice, if required. """
    if dont_deform_kpoint: return _kpoint
    if sum(abs(input-relaxed)) < 1e-11: return _kpoint
    return deform_kpoint(_kpoint, input, relaxed)
  def _line(start, end, density):
    """ Generator for creating points between two kpoints. """
    from numpy.linalg import norm

    distance = norm(_get_kpoint(end - start))
    nbkpt = int(max(1, density * distance - 1))
    stepsize = 1e0/float(nbkpt)
    _kpoints = [ float(i) * stepsize for i in range(1, nbkpt+1) ]
    for k in _kpoints: yield start + k * (end-start) 

  def _lines(endpoints, density):
    """ Generator for creating segments. """
    from numpy.linalg import norm
    assert len(endpoints) > 0, ValueError
    assert len(endpoints[0]) == 2, ValueError
    pos = 0
    yield pos, endpoints[0][0]
    for start, end in endpoints:
      last = start.copy()
      for _kpoint in _line(start, end, density):
        pos += norm(_get_kpoint(_kpoint-last))
        last = _kpoint.copy()
        yield pos, _kpoint


  # splits local communicator.
  color = comm.rank % pools
  local_comm = comm.split(color)
  # then computes different kpoints.
  results = []
  for i, (x, kpoint) in enumerate(_lines(kpoints, density)):
    # separates jobs into pools.
    if i % pools != color: continue

    # sets directory.
    directory = join(join(outdir, "band_structure"), "%i-%s" % (i, kpoint))
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

