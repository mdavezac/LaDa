#! /usr/bin/python
""" Module with different methods for vasp 

    Most of the generators in this module will not recompute the result if an
    appropriate output directory is found. In other words, this function can be
    called once to perform calculation, and then to perform output extractions.

    >>> # First create method.
    >>> some_method = SomeMethod(vasp)
    >>> # Then call the method.
    >>> some_method(structure)
    >>> # do something else
    >>>
    >>> # Then  analyze each step.
    >>> for extract in some_method.generator(structure):
    >>>   print "Where: ", extract.directory
    >>>   print extract.total_energy

    The output extraction object will be the output the vasp callable.
"""
from .extract import ExtractDFT, MassExtract
__docformat__ = "restructuredtext en"

class RelaxExtract(ExtractDFT):
  """ Extractor class for vasp relaxations. """
  class IntermediateMassExtract(MassExtract):
    """ Focuses on intermediate steps. """
    def __iter_alljobs__(self):
      """ Goes through all directories with an OUTVAR. """
      from glob import iglob
      from os.path import relpath, join
      from itertools import chain

      for dir in chain( iglob(join(join(self.rootdir, 'relax_cellshape'), '*/')),
                        iglob(join(join(self.rootdir, 'relax_ions'), '*/'))):
        try: result = ExtractDFT(dir[:-1])
        except: continue
        yield join('/', relpath(dir[:-1], self.rootdir)), result
        
  @property
  def details(self):
    """ Intermediate steps. """
    if '_details' not in self.__dict__:
      from os.path import exists
      if not exists(self.directory): return None
      self.__dict__['_details'] = None
      self._details = self.IntermediateMassExtract(self.directory)
      """ List of intermediate calculation extractors. """
    return self._details
  
  def iterfiles(self, **kwargs):
    """ Iterates over input/output files. """
    from itertools import chain
    for file in chain( super(EpiExtract, self).iterfiles(**kwargs),
                       self.details.iterfiles(**kwargs) ): yield file


class RelaxCellShape(object):
  """ Functor for cell-shape relaxation.
  
      Since vasp is a plane-wave code, cell-relaxation are never quite accurate.
      This functional keeps working until convergence is achieved, and then
      creates a static calculation.
  """
  Extract = RelaxExtract
  """ Extractor for relaxation functional. """
  def __init__( self, vasp, relaxation="volume ionic cellshape",\
                first_trial=None, maxiter=0, keep_steps=True ):
    """ Initializes a cell-shape relaxation. 
        
        :Parameters:
          vasp : lada.vasp.Vasp
            The functional with which to perform cell relaxation.
          relaxation
            Degrees of freedom to relax. Can be:
            
            - volume
            - ionic
            - cellshape

            Note that using this method with only *ionic* relaxation is a waste
            of computer ressources.
          first_trial : dict or None
            Does nothing if *None*.  A dictionary of parameters to apply for
            the first run of the functional. It may be interesting to first get a
            qualitative relaxation with less computationally intensive
            parameters of the vasp functional. 
          maxiter : int
            Maximum number of iterations before bailing out. Zero or negative
            number means infinit iterations.
          keep_steps : bool
	    If False, directory ``relax_cellshape`` where intermediate results
            are kept will be deleted once static are finished (and successfull).
    """
    super(RelaxCellShape, self).__init__()
    self.vasp = vasp
    """ The functional with which to perform cell relaxation. see `__init__`. """
    self.relaxation = relaxation
    """ Degrees of freedom to relax. see `__init__`. """
    self.first_trial = first_trial
    """ Dictionary of parameters for the first run of the functional. see `__init__`. """
    self.maxiter = maxiter
    """ Maximum number of iterations before bailing out. """
    self.keep_steps = keep_steps
    """ Whether or not to keep intermediate results. """


  def generator(self, structure, outdir=None, comm=None, **kwargs ):
    """ Performs a vasp relaxation, yielding each result.
    
        The result from the `vasp` calculations are yielded at each step. This
        makes it easy to analyse results during or after the run.

        >>> relaxor = RelaxCellShape(vasp)
        >>> for output in relaxor.generator(structure): 
        >>>   print output.total_energy

        :Parameters:
          structure
            The structure to relax with `vasp`.
          outdir
            Output directory passed on to the `vasp` functional.
          comm : mpi.Communicator or None
            MPI communicator passed on to the `vasp` functional.
          kwargs 
            Other keywords will overide attributes of this instance of
            `RelaxCellShape` (though for this run only. This function is
            stateless) if they are named after attributes of `RelaxCellShape`.
            Otherwise, the keywords are passed on to the `vasp` functional.
    """
    from copy import deepcopy
    from os import getcwd
    from os.path import join
    from shutil import rmtree
    from ..opt import RelativeDirectory
    from ..crystal import vasp_ordered
    from ..mpi import Communicator, world

    # make this function stateless.
    vasp = deepcopy(kwargs.pop("vasp", self.vasp))
    structure = vasp_ordered(structure, attributes=["magmom"])
    str_dict = deepcopy(structure.__dict__)
    first_trial = kwargs.pop("first_trial", self.first_trial)
    maxiter = kwargs.pop("maxiter", self.maxiter)
    keep_steps = kwargs.pop("keep_steps", self.keep_steps)
    outdir = getcwd() if outdir is None else RelativeDirectory(outdir).path
    nofail = kwargs.pop('nofail', False)
    if first_trial is None: first_trial = {}

    # convergence criteria and behavior.
    convergence = kwargs.get('convergence', getattr(self, 'convergence', self.vasp.ediffg))
    if convergence is None: convergence = 1e1 * self.vasp.ediff * float(len(structure.atoms))
    elif hasattr(convergence, "__call__"): pass
    elif convergence > 0: convergence *= float(len(structure.atoms))
    if convergence > 0 and convergence < self.vasp.ediff: 
      raise ValueError("Energy convergence criteria ediffg({0}) is smaller than ediff({1})."\
                       .format(self.vasp.ediffg, self.vasp.ediff))
    if hasattr(convergence, "__call__"):
      def is_converged(extractor):  
        if extractor is None: return True
        if not extractor.success: raise RuntimeError("VASP calculation did not succeed.")
        i = int(extractor.directory.split('/')[-1]) + 1
        if getattr(self, 'minrelsteps', i) > i: return False 
        return convergence(extractor)
    else:
      if convergence > 0e0:
        def is_converged(extractor):
          if extractor is None: return True
          if not extractor.success: raise RuntimeError("VASP calculation did not succeed.")
          i = int(extractor.directory.split('/')[-1]) + 1
          if getattr(self, 'minrelsteps', i) > i: return False 
          if extractor.total_energies.shape[0] < 2: return True
          return abs(extractor.total_energies[-2] - extractor.total_energies[-1:]) < convergence
      else:
        def is_converged(extractor):
          from numpy import max, abs, all
          if extractor is None: return True
          if not extractor.success: raise RuntimeError("VASP calculation did not succeed.")
          i = int(extractor.directory.split('/')[-1]) + 1
          if getattr(self, 'minrelsteps', i) > i: return False
          return all(max(abs(output.forces)) < abs(convergence))


    comm = Communicator(comm if comm is not None else world)

    # updates vasp as much as possible.
    if "set_relaxation" in kwargs: 
      from warnings import warn
      warn( DeprecationWarning("set_relaxation is deprecated. Please use relaxation."),\
            stacklevel=2 )
    vasp.relaxation = kwargs.pop("relaxation", kwargs.pop("set_relaxation", self.relaxation))
    for key in kwargs.keys():
      if hasattr(vasp, key): setattr(vasp, key, kwargs.pop(key))
     
    # does not run code. Just creates directory.
    if kwargs.pop("norun", False): 
      this = RelaxCellShape(vasp, vasp.relaxation, first_trial, maxiter)
      yield this._norun(structure, outdir=outdir, comm=comm, **kwargs)
      return

    # number of restarts.
    nb_steps, output = 0, None
   
    # sets parameter dictionary for first trial.
    if first_trial is not None:
      params = kwargs.copy()
      params.update(first_trial)
    else: params = kwargs
    comm.barrier()

    
    # performs relaxation calculations.
    while (maxiter <= 0 or nb_steps < maxiter) and vasp.relaxation.find("cellshape") != -1:
      # performs initial calculation.   
      output = vasp\
               (\
                 structure,
                 outdir = join(outdir, join("relax_cellshape", str(nb_steps))),
                 comm=comm,
                 restart = output,
                 **params
               )
      structure = output.structure
      structure.__dict__.update(str_dict)
      yield output
      assert output.success, RuntimeError("VASP calculations did not complete.")
      
      nb_steps += 1
      if nb_steps == 1 and len(first_trial) != 0: params = kwargs; continue
      # check for convergence.
      if is_converged(output): break;

    # Does not perform ionic calculation if convergence not reached.
    if nofail == False and is_converged(output) == False: 
      raise RuntimeError("Could not converge cell-shape in {0} iterations.".format(maxiter))

    # performs ionic calculation. 
    while (maxiter <= 0 or nb_steps < maxiter + 1) and vasp.relaxation.find("ionic") != -1:
      output = vasp\
               (\
                 structure, 
                 outdir = join(outdir, join("relax_ions", str(nb_steps))),
                 comm=comm,
                 relaxation = "ionic",
                 restart = output,
                 **params
               )
      structure = output.structure
      structure.__dict__.update(str_dict)
      yield output
      assert output.success, RuntimeError("VASP run did not succeed.")

      nb_steps += 1
      if nb_steps == 1 and len(first_trial) != 0: params = kwargs; continue
      # check for convergence.
      if is_converged(output): break;

    # Does not perform static calculation if convergence not reached.
    if nofail == False and is_converged(output) == False: 
      raise RuntimeError("Could not converge ions in {0} iterations.".format(maxiter))

    # performs final calculation outside relaxation directory. 
    output = vasp\
             (\
               structure, \
               outdir = outdir,\
               comm=comm,\
               relaxation = "static",\
               restart = output, \
               **kwargs\
             )
    yield output

    if output.success and (not keep_steps) and comm.is_root:
      rmtree(join(outdir, "relax_cellshape"))
      rmtree(join(outdir, "relax_ions"))

  def __call__(self, structure, outdir=None, comm=None, overwrite=False, **kwargs):
    """ Performs a vasp relaxation. 

        :Parameters:
          structure : lada.crystal.Structure
            The structure to relax with `vasp`.
          outdir
            Output directory passed on to the `vasp` functional.
          comm : mpi.communicator or None
            MPI communicator passed on to the `vasp` functional.
          overwrite : bool
	    Wether to perform the calculation no matter what, or whether to
            check if results already exist.
          kwargs 
            Other keywords will overide attributes of this instance of
            `RelaxCellShape` (though for this run only. This function is
            stateless) if they are named after attributes of `RelaxCellShape`.
            Otherwise, the keywords are passed on to the `vasp` functional.

        The plane wave basis depends upon cell-shape. Hence, the basis goes out
        of sync with the actual structure during cell-shape relaxation.
        Convergence can only be achieved by restarting the calculation. And
        eventually performing a static calculation.

        If you want to examine the result of each and every vasp calculation,
        use `generator` instead.
    """ 
    from os import getcwd
    from ..opt import RelativeDirectory
    from ..mpi import Communicator, world

    outdir = getcwd() if outdir is None else RelativeDirectory(outdir).path
    comm = Communicator(comm if comm is not None else world)
    if not overwrite:
      extract = self.Extract(outdir, comm=comm)
      if extract.success: return extract

    for output in self.generator(structure, outdir=outdir,
                                 comm=comm, overwrite=overwrite, **kwargs): pass
    return output
  
  def _norun(self, *args, **kwargs):
    """ Just creates directory for debugging. """
    from ..opt.changedir import Changedir
    from ..mpi import Communicator

    comm = Communicator(kwargs.get('comm', None))
    if comm.is_root:
      # creates a file describing the relaxation parameters.
      with Changedir(kwargs["outdir"]) as pwd:
        with open("relax_cell_shape_parameters", "w") as file:
          file.write( "self.relaxation = %s\n" % (repr(self.relaxation)) )
          file.write( "self.maxiter = %s\n" % (repr(self.maxiter)) )
          file.write( "self.first_trial = %s\n" % (repr(self.first_trial)) )

    # Now have vasp do a fake run to create anything it does create.
    kwargs["norun"] = True
    return self.vasp(*args, **kwargs) 

  def __repr__(self):
    """ Returns a python script describing this instance. """
    from re import sub
    string = "# VASP functional.\n"
    string += sub("(?<!\.)functional", "vasp_functional", repr(self.vasp))
    result = "from %s import %s\n\n" % (self.__class__.__module__, self.__class__.__name__)
    result += "# VASP Functional to relax cell-shape, volume, etc.\n"
    result += "functional = %s(vasp_functional)\n" % (self.__class__.__name__)
    result += "functional.relaxation = %s\n" % (repr(self.relaxation))
    result += "functional.first_trial = %s\n" % (repr(self.first_trial))
    result += "functional.maxiter = %s\n\n" % (repr(self.maxiter))
    return string + "\n" + result

def epi_relaxation( vasp, structure, outdir=None, comm=None,\
                    direction=[0,0,1], epiconv = 1e-4, final=None,
                    **kwargs ):
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
      :param dict final:
        parameters to change for final static calculation.
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
    for a in newstruct.atoms: # keep fractional coordinates constant.
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
  stress_direction = 1.0 if component(allcalcs[-1].stress) > 0e0 else -1.0
  xend = 0.1 if stress_direction > 0e0 else -0.1
  # compute xend value.
  eend = function(xend)
  # make sure xend is on other side of stress tensor sign.
  while stress_direction * component( allcalcs[-1].stress ) > 0e0:
    xstart, estart = xend, eend
    xend += 0.1 if stress_direction > 0e0 else -0.1
    eend = function(xend)
  
  # now we have a bracket. We start bisecting it.
  def convergence(a, b):
    """ Convergence criteria. """
    from numpy import dot
    from numpy.linalg import inv
    if epiconv > 0e0:
      return abs(a.total_energy - b.total_energy) \
             > epiconv * float(len(structure.atoms))
    else:
      avec = dot(a.structure.cell.T, direction) * a.structure.scale
      bvec = dot(b.structure.cell.T, direction) * b.structure.scale
      return any(abs(avec-bvec) > abs(epiconv))
      
     
  while convergence(estart, eend):
    xmid = 0.5 * (xend + xstart)
    emid = function(xmid)
    if stress_direction * component(emid.stress) > 0e0: xstart, estart = xmid, emid
    else: xend, eend = xmid, emid

  # last two calculation: relax mid-point of xstart, xend, then  perform static.
  xmid = 0.5 * (xend + xstart)
  emid = function(xmid)
  args = kwargs.copy()
  if final is not None: args.update(final)
  result = vasp( change_structure(xmid),
                 relaxation = "static",
                 outdir = outdir,
                 restart = allcalcs[-1],
                 comm = comm, **args )
  return result

# def epi_relaxation( vasp, structure, outdir=None, comm=None,\
# 		    direction=[0,0,1], epiconv = 1e-4, final=None,
# 		    **kwargs ):
#   """ Performs epitaxial relaxation in given direction. 
#   
#       Performs a relaxation for an epitaxial structure on a virtual substrate.
#       The external (cell) coordinates of the structure can only relax in the
#       growth/epitaxial direction. Internal coordinates (ions), however, are
#       allowed to relax in whatever direction. 
#       
#       Since VASP does not intrinsically allow for such a relaxation, it is
#       performed by chaining different vasp calculations together. The
#       minimization procedure itself is the secant method, enhanced by the
#       knowledge of the stress tensor. The last calculation is static, for
#       maximum accuracy.
# 
#       :param vasp: 
# 	:py:class:`Vasp <lada.Vasp>` functional with wich to perform the
# 	relaxation.
#       :param structure:
# 	:py:class:`Structure <lada.crystal.Structure>` for which to perform the
# 	relaxation.
#       :param str outdir: 
# 	Directory where to perform calculations. If None, defaults to current
# 	working directory. The intermediate calculations are stored in the
# 	relax_ions subdirectory.
#       :param comm: 
# 	Communicator with which to perform actual vasp calls. 
#       :param direction:
# 	Epitaxial direction. Defaults to [0, 0, 1].
#       :param float epiconv: 
# 	Convergence criteria of the total energy.
#       :param dict final:
# 	parameters to change for final static calculation.
#   """
#   from os import getcwd
#   from os.path import join
#   from copy import deepcopy
#   from numpy.linalg import norm
#   from numpy import array, dot
#   from lada.vasp.incar import PartialRestart
# 
#   direction = array(direction, dtype='float64') / norm(direction)
#   if outdir is None: outdir = getcwd()
# 
#   # creates relaxation functional.
#   vasp = deepcopy(vasp)
#   kwargs.pop('relaxation', None)
#   vasp.relaxation = 'ionic'
#   vasp.encut = 1.4
#   if 'encut' in kwargs: vasp.encut = kwargs.pop('encut')
#   if 'ediff' in kwargs: vasp.ediff = kwargs.pop('ediff')
#   if vasp.ediff < epiconv: vasp.ediff = epiconv * 1e-2
#   vasp.restart = PartialRestart(None)
#   kwargs['isif'] = 2
# 
#   allcalcs = []
#   def change_structure(rate):
#     """ Creates new structure with input change in c. """
#     from numpy.linalg import inv
#     if len(allcalcs) != 0: orig = allcalcs[-1].structure
#     else: orig = structure
#     newstruct = orig.copy()
#     cell = structure.cell
#     for i in xrange(3):
#       cell[:, i] += dot(structure.cell[:, i], direction) * rate * direction
#     newstruct.cell = cell
#     for a in newstruct.atoms: # keep fractional coordinates constant.
#       a.pos = dot(cell, dot(inv(orig.cell), a.pos))
#     return newstruct
# 
#   def component(stress):
#     """ Returns relevant stress component. """
#     return dot(dot(direction, stress), direction)
# 
#   def function(x):
#     """ Computes total energy for input change in c direction. """
#     e = vasp( change_structure(x),
# 	      outdir = join(outdir, "relax_ions/{0:0<12.10}".format(x)),
# 	      comm = comm,
# 	      restart = None if len(allcalcs) == 0 else allcalcs[-1],
# 	      **kwargs )
#     if not e.success:
#       raise RuntimeError("Vasp calculation in {0} did not complete.".format(e.directory))
#     allcalcs.append(e)
#     return e
# 
#   # Tries and find a bracket for minimum. 
#   # To do this, we start from current structure, look at stress in relevant
#   # direction for the direction in which to search, and expand/contract in that direction.
#   xstart = 0.0
#   estart = function(xstart)
#   sstart = component(estart.stress)
#   # then checks stress for actual direction to look at.
#   stress_direction = 1.0 if component(allcalcs[-1].stress) > 0e0 else -1.0
#   xend = 0.1 if stress_direction > 0e0 else -0.1
#   # compute xend value.
#   eend = function(xend)
#   send = component(eend.stress)
# 
#   # Now minimizes using gradient
#   if abs(send) > abs(sstart):  emid, xmid, smid = estart, xstart, sstart
#   else:  emid, xmid, smid = eend, xend, send
#   while abs(smid) > epiconv * len(structure.atoms):
#     if abs(xend-xstart) < 1e-8: break
#     send, sstart = component(eend.stress), component(estart.stress)
#     xmid = xend - send / (send - sstart) * (xend - xstart)
#     emid = function(xmid)
#     smid = component(emid.stress)
#     if abs(send) > abs(sstart): eend, xend, send = emid, xmid, smid
#     else: estart, xstart, sstart = emid, xmid, smid
#   
#   # last calculation: static.
#   if abs(send) > abs(sstart): xend, send = xstart, sstart
#   if abs(smid) > abs(send): xend = xmid
#   args = kwargs.copy()
#   if final is not None: args.update(final)
#   return vasp( change_structure(xend),
# 	       relaxation = "static",
# 	       outdir = outdir,
# 	       restart = allcalcs[-1],
# 	       comm = comm, **args )
epi_relaxation.Extract = RelaxExtract
