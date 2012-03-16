""" Relaxation Methods

    An accurately strain-relaxed VASP calculation requires multiple restarts.
    The reasons for this lies in that the plane-wave basis is determined at the
    start of any particular VASP run. Hence, the basis is incorrect if the
    cell-shape changes during the run. The same can be said of real-space
    pseudo-potential grids when relaxing ionic positions. 

    This module contains methods to chain together VASP calculations until a
    fully relaxed structure is obtained.
"""
__docformat__ = "restructuredtext en"
__all__ = ['relax', 'iter_relax', 'epitaxial', 'iter_epitaxial']
from .extract import Extract

def iter_relax( vasp, structure, outdir=None, first_trial=None,
                maxcalls=10, keepsteps=True, nofail=False, 
                convergence=None, relaxation = "volume cellshape ionic", 
                minrelsteps=-1, **kwargs ):
  """ Iterator over calls to VASP during relaxation.
  
      This generator iterates over successive VASP calculations until a fully
      relaxed structure is obtained. Its last calculation is *static*, ensuring
      that the final electronic structure accurately represents the relaxed
      structure.

      The full process is to first relax the cell-shape (and internal degrees of
      freedom upon request) until convergence is achieved, as determined by the
      difference in total energies (see the keyword argument ``convergence``)
      within the current VASP run. Subsequent runs keep the cell-shape constant
      while allowing ionic degrees of freedom to relax, until the same
      convergence criteria is achieved. Finally, a static calculation is
      performed.

      It is possible to bypass cell-shape relaxations and perform only
      ionic-relaxations.

      :param vasp:
        :py:class:`Vasp <lada.vasp.functional.Vasp>` object with which to
        perform the relaxation.
      :param structure:
        :py:class:`Structure <lada.crystal.Structure>` object for which to
        perform the relaxation.
      :param outdir:
        Directory where to perform the calculations. Defaults to current
        working directory. The actual calculations are stored within the
        *relax* subdirectory.
      :param dict first_trial:
        Holds parameters which are used only for the very first VASP
        calculation. It can be used to accelerate the first step of the
        relaxation if starting far from the optimum. Defaults to empty
        dictionary.
      :param int maxcalls:
        Maximum number of calls to VASP before aborting. Defaults to 10.
      :param bool keepsteps:
        If true, intermediate steps are kept. If False, intermediate steps are
        erased.
      :param string relaxation:
        Degrees of freedom to relax. It should be either "cellshape" or "ionic"
        or both.
      :param bool nofail:
        If True, will not fail if convergence is not achieved. Just keeps going. 
        Defaults to False.
      :param convergence:
        Convergence criteria. If ``minrelsteps`` is positive, it is only
        checked after ``minrelsteps`` have been performed.

        * None: defaults to ``vasp.ediff * 1e1``
        * positive real number: energy convergence criteria in eV. 
        * negative real number: force convergence criteria in eV/angstrom. 
        * callable: Takes an extraction object as input. Should return True if
          convergence is achieved and False otherwise.
      :param int minrelsteps:
        Fine tunes how convergence criteria is applied.
        
        * positive: at least ``minrelsteps`` calls to VASP are performed before
          checking for convergence. If ``relaxation`` contains "cellshape",
          then these calls occur during cellshape relaxation. If it does not,
          then the calls occur during the ionic relaxations. The calls do count
          towards ``maxcalls``.
        * negative (default): argument is ignored.
      :param kwargs:
        Other parameters are applied to the input :py:class:`Vasp
        <lada.vasp.functional.Vasp>` object.

      :return: At each step, yields an extraction object if the relevant VASP
               calculation already exists. Otherwise, it yields a
               :py:class:`Program <lada.misc.Program>` object detailing the
               call to the external VASP program.

      .. seealso:: :py:func:`execute_program <lada.misc.execute_program>`
  """
  from re import sub
  from copy import deepcopy
  from os import getcwd
  from os.path import join
  from shutil import rmtree
  from ..misc import RelativePath
  from ..crystal import vasp_ordered

  # make this function stateless.
  vasp = deepcopy(vasp)
  relaxed_structure = structure.copy()
  if first_trial is None: first_trial = {}
  outdir = getcwd() if outdir is None else RelativePath(outdir).path

  # convergence criteria and behavior.
  is_converged = _get_is_converged( vasp, relaxed_structure, convergence=convergence,
                                    minrelsteps=minrelsteps, **kwargs)

  # number of restarts.
  nb_steps, output = 0, None
 
  # sets parameter dictionary for first trial.
  if first_trial is not None:
    params = kwargs.copy()
    params.update(first_trial)
  else: params = kwargs
  
  # performs relaxation calculations.
  while (maxcalls <= 0 or nb_steps < maxcalls) and relaxation.find("cellshape") != -1:
    # performs initial calculation.   
    for u in vasp.iter\
             (\
               relaxed_structure,
               outdir = join(outdir, join("relax_cellshape", str(nb_steps))),
               restart = output,
               relaxation = relaxation,
               **params
             ): yield u

    output = vasp.Extract(join(outdir, join("relax_cellshape", str(nb_steps))))
    assert output.success, RuntimeError("VASP calculations did not complete.")
    relaxed_structure = output.structure
    
    nb_steps += 1
    if nb_steps == 1 and len(first_trial) != 0: params = kwargs; continue
    # check for convergence.
    if is_converged(output): break;

  # Does not perform ionic calculation if convergence not reached.
  if nofail == False and is_converged(output) == False: 
    raise RuntimeError("Could not converge cell-shape in {0} iterations.".format(maxcalls))

  # performs ionic calculation. 
  while (maxcalls <= 0 or nb_steps < maxcalls + 1) and relaxation.find("ionic") != -1:
    for u in vasp.iter\
             (\
               relaxed_structure, 
               outdir = join(outdir, join("relax_ions", str(nb_steps))),
               relaxation = "ionic",
               restart = output,
               **kwargs
             ): yield u

    output = vasp.Extract(join(outdir, join("relax_ions", str(nb_steps))))
    assert output.success, RuntimeError("VASP calculations did not complete.")
    relaxed_structure = output.structure

    nb_steps += 1
    if nb_steps == 1 and len(first_trial) != 0: params = kwargs; continue
    # check for convergence.
    if is_converged(output): break;

  # Does not perform static calculation if convergence not reached.
  if nofail == False and is_converged(output) == False: 
    raise RuntimeError("Could not converge ions in {0} iterations.".format(maxcalls))

  # performs final calculation outside relaxation directory. 
  for u in vasp.iter\
           (\
             relaxed_structure, \
             outdir = outdir,\
             relaxation = "static",\
             restart = output, \
             **kwargs\
           ): yield u
  output = vasp.Extract(outdir)
  assert output.success, RuntimeError("VASP calculations did not complete.")

  # replace initial structure with that with which this function was called.
  with output.__outcar__() as file:
    filename = file.name
    string = sub(  '#+ INITIAL STRUCTURE #+\n((.|\n)*)\n#+ END INITIAL STRUCTURE #+',
                   """################ INITIAL STRUCTURE ################\n"""\
                   """from {0.__class__.__module__} import {0.__class__.__name__}\n"""\
                   """structure = {1}\n"""\
                   """################ END INITIAL STRUCTURE ################\n"""\
                   .format(structure, repr(structure).replace('\n', '\n            ')),
                   file.read() )
  with open(filename, 'w') as file: file.write(string)

  if output.success and (not keepsteps):
    rmtree(join(outdir, "cellshape"))
    rmtree(join(outdir, "ions"))
iter_relax.Extract = Extract
""" Extraction method for relaxation runs. """

def relax( vasp, structure, outdir=None, first_trial=None,
           maxcalls=10, keepsteps=True, nofail=False, 
           convergence=None, relaxation = "volume cellshape ionic", 
           minrelsteps=-1, **kwargs ):
  """ Performs relaxation of an input structure using :py:class:`Vasp`.  """
  from os.path import join
  from ..misc import execute_program
  result = None
  for program in iter_relax( vasp, structure, outdir, first_trial=first_trial,
                             maxcalls=maxcalls, keepsteps=keepsteps, nofail=nofail,
                             relaxation=relaxation, **kwargs ):
    # iterator may yield the result from a prior successfull run. 
    if getattr(program, 'success', False):
      result = program
      continue
    # otherwise, it should yield a Program tuple to execute.
    execute_program(program, append=False, **kwargs)
    result = vasp.Extract(outdir)
  return result

relax.__doc__ += iter_relax.__doc__[iter_relax.__doc__.find('\n'):\
                                    iter_relax.__doc__.find(':return')]\
                                   .replace('generator', 'method').replace('\n      ', '\n')\
                 + "\n:return: An extraction object pointing to the final static calculation.\n"\
                 + "\n.. seealso:: `iter_relax`\n\n"
relax.Extract = iter_relax.Extract
""" Extraction method for relaxation runs. """

def _get_is_converged(vasp, structure, convergence=None, minrelsteps=-1, **kwargs):
  """ Returns convergence function. """
  # tries and devine the convergence criteria from the input.
  if convergence is None: convergence = 1e1 * vasp.ediff
  elif hasattr(convergence, "__call__"): pass
  elif convergence > 0: convergence *= float(len(structure))
  if convergence > 0 and convergence < vasp.ediff: 
    raise ValueError("Energy convergence criteria ediffg({0}) is smaller than ediff({1})."\
                     .format(vasp.ediffg, vasp.ediff))
  # creates a convergence function.
  if hasattr(convergence, "__call__"):
    def is_converged(extractor):  
      if extractor is None: return True
      if not extractor.success: raise RuntimeError("VASP calculation did not succeed.")
      i = int(extractor.directory.split('/')[-1]) + 1
      if minrelsteps > 0 and minrelsteps > i: return False
      return convergence(extractor)
  elif convergence > 0e0:
    def is_converged(extractor):
      if extractor is None: return True
      if not extractor.success: raise RuntimeError("VASP calculation did not succeed.")
      i = int(extractor.directory.split('/')[-1]) + 1
      if minrelsteps > 0 and minrelsteps > i: return False
      if extractor.total_energies.shape[0] < 2: return True
      return abs(extractor.total_energies[-2] - extractor.total_energies[-1:]) < convergence
  else:
    def is_converged(extractor):
      from numpy import max, abs, all
      if extractor is None: return True
      if not extractor.success: raise RuntimeError("VASP calculation did not succeed.")
      i = int(extractor.directory.split('/')[-1]) + 1
      if minrelsteps > 0 and minrelsteps > i: return False
      return all(max(abs(output.forces)) < abs(convergence))
  return is_converged

def iter_epitaxial(vasp, structure, outdir=None, direction=[0,0,1], epiconv = 1e-4,
                   initstep=0.05, **kwargs):
  """ Performs epitaxial relaxation in given direction. 
  
      This generator iterates over successive VASP calculations until an
      epitaxially relaxed structure is obtained.  The external (cell)
      coordinates of the structure can only relax in the growth/epitaxial
      direction. Internal coordinates (ions), however, are allowed to relax in
      whatever direction. 
      
      Since VASP does not intrinsically allow for such a relaxation, it is
      performed by chaining different vasp calculations together. The
      minimization procedure itself is the secant method, enhanced by the
      knowledge of the stress tensor. The last calculation is static, for
      maximum accuracy.

      :param vasp: 
        :py:class:`Vasp <lada.vasp.functional.Vasp>` functional with wich to
        perform the relaxation.
      :param structure:
        :py:class:`Structure <lada.crystal.Structure>` for which to perform the
        relaxation.
      :param str outdir: 
        Directory where to perform calculations. If None, defaults to current
        working directory. The intermediate calculations are stored in the
        relax_ions subdirectory.
      :param direction:
        Epitaxial direction. Defaults to [0, 0, 1].
      :param float epiconv: 
        Convergence criteria of the total energy.
      
      :return: At each step, yields an extraction object if the relevant VASP
               calculation already exists. Otherwise, it yields a
               :py:class:`Program <lada.misc.Program>` object detailing the
               call to the external VASP program.

      .. seealso:: :py:func:`execute_program <lada.misc.execute_program>`
  """
  from os import getcwd
  from os.path import join
  from copy import deepcopy
  from re import sub
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
  if kwargs.get('relaxation', None) is not None:
    vasp.relaxation = kwargs['relaxation']
  kwargs['relaxation'] = 2

  allcalcs = []
  def change_structure(rate):
    """ Creates new structure with input change in c. """
    from numpy.linalg import inv
    from numpy import outer, dot
    newstruct = structure.copy()
    strain = outer(direction, direction) * x 
    result.cell += dot(strain, structure.cell)
    for atom in result: atom.pos += dot(strain, atom.pos)
    return newstruct

  def component(stress):
    """ Returns relevant stress component. """
    return dot(dot(direction, stress), direction)

  # Tries and find a bracket for minimum. 
  # To do this, we start from current structure, look at stress in relevant
  # direction for the direction in which to search, and expand/contract in that direction.
  xstart = 0.0
  for u in vasp.iter( change_structure(xstart),
                      outdir = join(outdir, join("relax_ions", "{0:0<12.10}".format(xstart))),
                      restart = None if len(allcalcs) == 0 else allcalcs[-1],
                      **kwargs ): yield u
  estart = vasp.Extract(join(join(outdir, 'relax_ions'), '{0:0<12.10}'.format(xstart)))
  allcalcs.append(estart)
  
  # then checks stress for actual direction to look at.
  stress_direction = 1.0 if component(allcalcs[-1].stress) > 0e0 else -1.0
  xend = initstep if stress_direction > 0e0 else -initstep
  # compute xend value.
  for u in vasp.iter( change_structure(xend),
                      outdir = join(outdir, join("relax_ions", "{0:0<12.10}".format(xend))),
                      restart = None if len(allcalcs) == 0 else allcalcs[-1],
                      **kwargs ): yield u
  eend = vasp.Extract(join(join(outdir, 'relax_ions'), '{0:0<12.10}'.format(xend)))
  allcalcs.append(eend)
  # make sure xend is on other side of stress tensor sign.
  while stress_direction * component( allcalcs[-1].stress ) > 0e0:
    xstart, estart = xend, eend
    xend += initstep if stress_direction > 0e0 else -initstep
    for u in vasp.iter( change_structure(xend),
                        outdir = join(outdir, join("relax_ions", "{0:0<12.10}".format(xend))),
                        restart = None if len(allcalcs) == 0 else allcalcs[-1],
                        **kwargs ): yield u
    eend = vasp.Extract(join(join(outdir, 'relax_ions'), '{0:0<12.10}'.format(xend)))
    allcalcs.append(eend)
  
  # now we have a bracket. We start bisecting it.
  while abs(estart.total_energy - eend.total_energy) > epiconv * float(len(structure)):
    xmid = 0.5 * (xend + xstart)
    for u in vasp.iter( change_structure(xmid),
                        outdir = join(outdir, join("relax_ions", "{0:0<12.10}".format(xmid))),
                        restart = None if len(allcalcs) == 0 else allcalcs[-1],
                        **kwargs ): yield u
    emid = vasp.Extract(join(join(outdir, 'relax_ions'), '{0:0<12.10}'.format(xmid)))
    allcalcs.append(emid)
    if stress_direction * component(emid.stress) > 0: xstart, estart = xmid, emid
    else: xend, eend = xmid, emid

  # last two calculation: relax mid-point of xstart, xend, then  perform static.
  efinal = eend if estart.total_energy > eend.total_energy else estart
  kwargs['relaxation'] = 'static'
  for u in vasp.iter(efinal.structure, outdir=outdir, restart=efinal, **kwargs): yield u
  final = vasp.Extract(outdir)

  # replace initial structure with that with which this function was called.
  with final.__outcar__() as file:
    filename = file.name
    string = sub(  '#+ INITIAL STRUCTURE #+\n((.|\n)*)\n#+ END INITIAL STRUCTURE #+',
                   """################ INITIAL STRUCTURE ################\n"""\
                   """from {0.__class__.__module__} import {0.__class__.__name__}\n"""\
                   """structure = {1}\n"""\
                   """################ END INITIAL STRUCTURE ################\n"""\
                   .format(structure, repr(structure).replace('\n', '\n            ')),
                   file.read() )
  with open(filename, 'w') as file: file.write(string)
iter_epitaxial.Extract = Extract
""" Extraction method for epitaxial relaxation runs. """

def epitaxial(vasp, structure, outdir=None, direction=[0,0,1], epiconv = 1e-4, 
              initstep=0.05, **kwargs):
  """ Iterates over calls to VASP during epitaxial relaxation. """
  from os.path import join
  from ..misc import execute_program
  result = None
  for program in iter_epitaxial( vasp, structure, outdir, direction=direction,
                                 epiconv=epiconv, **kwargs ): 
    # iterator may yield the result from a prior successfull run. 
    if getattr(program, 'success', False):
      result = program
      continue
    # otherwise, it should yield a Program tuple to execute.
    execute_program(program, append=False, **kwargs)
    result = vasp.Extract(outdir)

  return result
epitaxial.Extract = iter_epitaxial.Extract
""" Extraction method for epitaxial relaxation runs. """
epitaxial.__doc__\
    += iter_epitaxial.__doc__[iter_epitaxial.__doc__.find('\n'):\
                              iter_epitaxial.__doc__.find(':return')]\
                             .replace('generator', 'method').replace('\n      ', '\n')\
                 + "\n:return: An extraction object pointing to the final static calculation.\n"\
                 + "\n.. seealso:: `iter_epitaxial`\n\n"
