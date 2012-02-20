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
__docformat__ = "restructuredtext en"

class RelaxCellShape(object):
  """ Functor for cell-shape relaxation.
  
      Since vasp is a plane-wave code, cell-relaxation are never quite accurate.
      This functional keeps working until convergence is achieved, and then
      creates a static calculation.
  """
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

  @property
  def Extract(self):
    """ Extraction class. """
    return self.vasp.Extract

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
    while maxiter <= 0 or nb_steps < maxiter and vasp.relaxation.find("cellshape") != -1:
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
    while maxiter <= 0 or nb_steps < maxiter + 1 and vasp.relaxation.find("ionic") != -1:
      output = vasp\
               (\
                 structure, 
                 outdir = join(outdir, join("relax_ions", str(nb_steps))),
                 comm=comm,
                 relaxation = "ionic",
                 restart = output,
                 **kwargs
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
    string = "# VASP functional.\n"
    string += repr(self.vasp).replace("functional", "vasp_functional") 
    result = "from %s import %s\n\n" % (self.__class__.__module__, self.__class__.__name__)
    result += "# VASP Functional to relax cell-shape, volume, etc.\n"
    result += "functional = %s(vasp_functional)\n" % (self.__class__.__name__)
    result += "functional.relaxation = %s\n" % (repr(self.relaxation))
    result += "functional.first_trial = %s\n" % (repr(self.first_trial))
    result += "functional.maxiter = %s\n\n" % (repr(self.maxiter))
    return string + "\n" + result
