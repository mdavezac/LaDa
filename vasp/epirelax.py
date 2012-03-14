def iter_epirelax( vasp, structure, outdir=None, comm=None,\
                   direction=[0,0,1], epiconv = 1e-4, **kwargs ):
  """ Performs epitaxial relaxation in given direction. 
  
      Performs a relaxation for an epitaxial structure on a virtual substrate.
      The external (cell) coordinates of the structure can only relax in the
      growth/epitaxial direction. Internal coordinates (ions), however, are
      allowed to relax in whatever direction. 
      
      Since VASP does not intrinsically allow for such a relaxation, it is
      performed by chaining different vasp calculations together. The
      minimization procedure itself is the secant method, enhanced by the
      knowledge of the stress tensor. The last calculation is static, for
      maximum accuracy.

      :param vasp: 
        :py:class:`Vasp <lada.Vasp>` functional with wich to perform the
        relaxation.
      :param structure:
        :py:class:`Structure <lada.crystal.Structure>` for which to perform the
        relaxation.
      :param str outdir: 
        Directory where to perform calculations. If None, defaults to current
        working directory. The intermediate calculations are stored in the
        relax_ions subdirectory.
      :param comm: 
        Communicator with which to perform actual vasp calls. 
      :param direction:
        Epitaxial direction. Defaults to [0, 0, 1].
      :param float epiconv: 
        Convergence criteria of the total energy.
  """
  from os import getcwd
  from os.path import join
  from copy import deepcopy
  from numpy.linalg import norm
  from numpy import array, dot
  from lada.vasp.incar import PartialRestart

  direction = array(direction, dtype='float64') / norm(direction)
  if outdir is None: outdir = getcwd()

  # creates relaxation functional.
  vasp = deepcopy(vasp)
  kwargs.pop('relaxation', None)
  vasp.relaxation = 'ionic'
  vasp.encut = 1.4
  if 'encut' in kwargs: vasp.encut = kwargs.pop('encut')
  if 'ediff' in kwargs: vasp.ediff = kwargs.pop('ediff')
  if vasp.ediff < epiconv: vasp.ediff = epiconv * 1e-2
  vasp.restart = PartialRestart(None)
  kwargs['isif'] = 2

  allcalcs = []
  def change_structure(rate):
    """ Creates new structure with input change in c. """
    from numpy.linalg import inv
    if len(allcalcs) != 0: orig = allcalcs[-1].structure
    else: orig = structure
    newstruct = orig.copy()
    cell = structure.cell
    for i in xrange(3):
      cell[:, i] += dot(structure.cell[:, i], direction) * rate * direction
    newstruct.cell = cell
    for a in newstruct: # keep fractional coordinates constant.
      a.pos = dot(cell, dot(inv(orig.cell), a.pos))
    return newstruct

  def component(stress):
    """ Returns relevant stress component. """
    return dot(dot(direction, stress), direction)

  def function(x):
    """ Computes total energy for input change in c direction. """
    e = vasp( change_structure(x),
              outdir = join(outdir, "relax_ions/{0:0<12.10}".format(x)),
              comm = comm,
              restart = None if len(allcalcs) == 0 else allcalcs[-1],
              **kwargs )
    if not e.success:
      raise RuntimeError("Vasp calculation in {0} did not complete.".format(e.directory))
    allcalcs.append(e)
    return e

  # Tries and find a bracket for minimum. 
  # To do this, we start from current structure, look at stress in relevant
  # direction for the direction in which to search, and expand/contract in that direction.
  xstart = 0.0
  estart = function(xstart)
  # then checks stress for actual direction to look at.
  stress_direction = 1.0 if component(allcalcs[-1].stress) else -1.0
  xend = 0.1 if stress_direction else -0.1
  # compute xend value.
  eend = function(xend)
  # make sure xend is on other side of stress tensor sign.
  while stress_direction * component( allcalcs[-1].stress ) > 0e0:
    xstart, estart = xend, eend
    xend += 0.1
    eend = function(xend)
  
  # now we have a bracket. We start bisecting it.
  while abs(estart.total_energy - eend.total_energy) > epiconv * float(len(structure.atoms)):
    xmid = 0.5 * (xend + xstart)
    emid = function(xmid)
    if component(emid.stress) > 0: xstart, estart = xmid, emid
    else: xend, eend = xmid, emid

  # last two calculation: relax mid-point of xstart, xend, then  perform static.
  xmid = 0.5 * (xend + xstart)
  emid = function(xmid)
  result = vasp( change_structure(xmid),
                 relaxation = "static",
                 outdir = outdir,
                 restart = allcalcs[-1],
                 comm = comm, **kwargs )
  return result
