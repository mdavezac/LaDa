#! python
""" Bandstructure plotting tools """
from sys import exit

def band_structure(escan, structure, kpoints, density, outdir=None, comm=None,\
                   do_relax=None, pools = 1, **kwargs):
  """ Returns eigenvalues for plotting bandstructure. """
  from os import getcwd
  from os.path import join, expanduser, abspath, exists
  from shutil import copyfile
  from boost.mpi import world, all_gather
  from numpy.linalg import norm
  from numpy import abs, sum
  from ..crystal import deform_kpoint
  from ..opt.changedir import Changedir

  # check/correct input arguments
  assert "do_genpot" not in kwargs,\
         ValueError("\"do_genpot\" is not an admissible argument of band_structure.")
  assert "do_escan" not in kwargs,\
         ValueError("\"do_escan\" is not an admissible argument of band_structure.")
  outdir = abspath(expanduser(outdir)) if outdir != None else getcwd()
  if do_relax == None: do_relax = escan.do_relax 
  if comm == None: comm = world
  if pools > comm.size: pools = comm.size

  # first computes vff and genpot.
  potdir = join(outdir, "band_structure")
  vffout = escan( structure, outdir=potdir, do_escan=False, do_genpot=True,\
                  do_relax=do_relax, comm = comm, **kwargs )
  
  # two functions required to continue.
  def _get_kpoint(_kpoint):
    """ Deforms kpoint to new lattice, if required. """
    if vffout.escan._dont_deform_kpoint: return _kpoint
    input, relaxed = structure.cell, vffout.structure.cell
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
    # copies POSCAR and POTCAR for reuse.
    with Changedir(directory, comm = local_comm) as cwd:
      POSCAR = escan._POSCAR + "." + str(local_comm.rank)
      if exists(join(potdir, POSCAR)): copyfile(join(potdir, POSCAR), POSCAR)
      POTCAR = escan._POTCAR + "." + str(local_comm.rank)
      if exists(join(potdir, POTCAR)): copyfile(join(potdir, POTCAR), POTCAR)
      cout = escan.vff._cout(comm)
      if exists(join(potdir, cout)): copyfile(join(potdir, cout), cout)

    # actually computes stuff.
    out = escan( structure, outdir=directory, kpoint=kpoint, do_relax=False,\
                 do_genpot=False, do_escan=True, comm = local_comm, **kwargs )
    # saves stuff
    eigenvalues = out.eigenvalues.copy()
    eigenvalues.sort()
    results.append( (x, kpoint, eigenvalues) )

  if comm.size > 1: # gathers and orders results.
    head_comm = comm.split(0 if local_comm.rank == 0 else 1)
    if local_comm.rank == 0:
      results = all_gather(head_comm, results)
      results = sorted((j for i in results for j in i), lambda a,b: a[0] < b[0])
      broadcast(local_comm, results, 0)
    else: results = broadcast(local_comm, None, 0) 
  return results

