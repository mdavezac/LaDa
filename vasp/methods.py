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
    self.Extract = self.vasp.Extract
    """ Extraction class. """
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
          comm : boost.mpi.communicator or None
            MPI communicator passed on to the `vasp` functional.
          kwargs 
            Other keywords will overide attributes of this instance of
            `RelaxCellShape` (though for this run only. This function is
            stateless) if they are named after attributes of `RelaxCellShape`.
            Otherwise, the keywords are passed on to the `vasp` functional.
    """
    from warnings import warn
    from copy import deepcopy
    from math import fabs 
    from os import getcwd
    from os.path import join
    from shutil import rmtree
    from ..opt import RelativeDirectory


    # make this function stateless.
    vasp = deepcopy(kwargs.pop("vasp", self.vasp))
    structure = deepcopy(structure)
    str_dict = deepcopy(structure.__dict__)
    first_trial = kwargs.pop("first_trial", self.first_trial)
    maxiter = kwargs.pop("maxiter", self.maxiter)
    keep_steps = kwargs.pop("keep_steps", self.keep_steps)
    outdir = getcwd() if outdir == None else RelativeDirectory(outdir).path
    ediffg = self.vasp.ediffg
    if ediffg == None: ediffg = 1e1 * self.vasp.ediff
    elif ediffg < self.vasp.ediff: 
      raise ValueError( "Parameter ediffg({0}) is smaller than ediffg({1}."\
                        .format(self.ediffg, self.vasp.ediff))
    ediffg *= 1.2 * float(len(structure.atoms))

    # updates vasp as much as possible.
    if "set_relaxation" in kwargs: 
      warn("set_relaxation is deprecated. Please use relaxation.", DeprecationWarning)
    vasp.relaxation = kwargs.pop("relaxation", kwargs.pop("set_relaxation", self.relaxation))
    for key in kwargs.keys():
      if hasattr(vasp, key): setattr(vasp, key, kwargs.pop(key))
     
    # does not run code. Just creates directory.
    if kwargs.pop("norun", False): 
      this = RelaxCellShape(vasp, relaxation, first_trial, maxiter)
      yield this._norun(structure, outdir=outdir, comm=comm, **kwargs)
      return

    # number of restarts.
    nb_steps, output = 0, None
   
    # sets parameter dictionary for first trial.
    if first_trial != None:
      params = kwargs.copy()
      params.update(first_trial)
    else: params = kwargs
    if comm != None: comm.barrier()
    
    # performs relaxation calculations.
    while maxiter <= 0 or nb_steps < maxiter and vasp.relaxation.find("cellshape") != -1:
      # performs initial calculation.   
      output = vasp\
               (\
                 structure,
                 outdir = join(outdir, join("relax_cellshape", str(nb_steps))),
                 comm=comm,
                 restart = output if nb_steps > 1 else None,
                 **params
               )
      structure = output.structure
      structure.__dict__.update(str_dict)
      yield output
      assert output.success, RuntimeError("VASP calculations did not complete.")
      
      nb_steps += 1
      if nb_steps == 1: params = kwargs; continue
      if output.total_energies.shape[0] < 2: break
      energies = output.total_energies[-2] - output.total_energies[-1:]
      if abs(energies) < ediffg: break

    # Does not perform ionic calculation if convergence not reached.
    if output != None:
      assert output.success, RuntimeError("VASP calculations did not complete.")
      if output.total_energies.shape[0] >= 2:
        energies = output.total_energies[-2] - output.total_energies[-1:]
        assert abs(energies) < ediffg, \
               RuntimeError("Could not converge cell-shape in {0} iterations.".format(maxiter))

    # performs ionic calculation. 
    while maxiter <= 0 or nb_steps < maxiter + 1 and vasp.relaxation.find("ionic") != -1:
      output = vasp\
               (\
                 structure, 
                 outdir = join(outdir, join("relax_ions", str(nb_steps))),
                 comm=comm,
                 relaxation = "ionic",
                 restart = output if nb_steps > 1 else None,
                 **kwargs
               )
      structure = output.structure
      structure.__dict__.update(str_dict)
      yield output
      assert output.success, RuntimeError("VASP run did not succeed.")

      nb_steps += 1
      if nb_steps == 1: params = kwargs; continue
      if output.total_energies.shape[0] < 2: break
      energies = output.total_energies[-2] - output.total_energies[-1:]
      if abs(energies) < ediffg: break

    # Does not perform static calculation if convergence not reached.
    if output != None:
      assert output.success, RuntimeError("VASP calculations did not complete.")
      if output.total_energies.shape[0] >= 2:
        energies = output.total_energies[-2] - output.total_energies[-1:]
        assert abs(energies) < ediffg, \
               RuntimeError("Could not converge ions in {0} iterations.".format(maxiter))

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

    if output.success and (not keep_steps):
      rmtree(join(outdir, "relax_cellshape"))
      rmtree(join(outdir, "relax_ions"))

  def __call__(self, structure, outdir=None, comm=None, overwrite=False, **kwargs):
    """ Performs a vasp relaxation. 

        :Parameters:
          structure : lada.crystal.Structure
            The structure to relax with `vasp`.
          outdir
            Output directory passed on to the `vasp` functional.
          comm : boost.mpi.communicator or None
            MPI communicator passed on to the `vasp` functional.
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
    from shutil import rmtree
    from os import getcwd
    from os.path import exists
    from ..opt import RelativeDirectory

    outdir = getcwd() if outdir == None else RelativeDirectory(outdir).path
    if not overwrite:
      extract = self.Extract(outdir, comm=None)
      if extract.success: return extract
    elif is_root and exists(outdir): rmtree(outdir)
    if comm != None: comm.barrier() # makes sure directory is not created by other proc!

    for output in self.generator(structure, outdir=outdir,
                                 comm=comm, overwrite=overwrite, **kwargs): pass
    return output
  
  def _norun(self, *args, **kwargs):
    """ Just creates directory for debugging. """
    from ..opt.changedir import Changedir

    is_root = True
    if "comm" in kwargs:
      comm = kwargs["comm"]
      is_root = True if comm == None else comm.rank == 0
    if is_root:
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
