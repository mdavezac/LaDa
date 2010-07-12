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
  def __init__(self, vasp, relaxation="volume ionic cellshape", first_trial=None, maxiter=0):
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

  def generator(self, structure, outdir=None, comm=None, **kwargs ):
    """ Performs a vasp relaxation, yielding each result.
    
        The result from the `vasp` calculations are yielded at each step. This
        makes it easy to analyse results during or after the run.

        >>> relaxor = RelaxCellShape(vasp)
        >>> for output in relaxor.generator(structure): 
        >>>   print output.energy

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
    from copy import deepcopy
    from math import fabs 
    from os import getcwd
    from os.path import join, exists


    # make this function stateless.
    vasp = kwargs.pop("vasp", self.vasp)
    structure = deepcopy(structure)
    set_relaxation = kwargs.pop("relaxation", self.relaxation)
    set_relaxation = kwargs.pop("set_relaxation", set_relaxation)
    first_trial = kwargs.pop("first_trial", self.first_trial)
    maxiter = kwargs.pop("maxiter", self.maxiter)
    outdir = kwargs.pop("outdir", getcwd())


    # number of restarts.
    nb_steps, olde = 0, None
   
    # sets parameter dictionary for first trial.
    if first_trial != None and nb_steps == 0:
      params = kwargs.copy()
      params.update(first_trial)
    else: param = kwargs
    
    # performs relaxation calculations.
    while maxiter <= 0 or nbsteps < maxiter:
      # performs initial calculation.
      output = vasp\
               (\
                 structure, \
                 outdir = join(outdir, "relax_cellshape_" + str(nb_steps)),\
                 comm=comm,\
                 set_relaxation = set_relaxation,\
                 **params
               )
      yield output
      
      nb_steps += 1
      params = kwargs
      olde = output.energy
      if nb_steps == 1: continue
      if fabs(output.energy - olde) < float(len(structure.atoms)) * self.vasp.ediff: break

    # Does not perform static calculation if convergence not reached.
    if fabs(output.energy - olde) > float(len(structure.atoms)) * self.vasp.ediff: 
      yield output 
    # performs final calculation outside relaxation directory. 
    output = vasp\
             (\
               structure, \
               outdir = outdir,\
               comm=comm,\
               set_relaxation = "static",\
               **kwargs\
             )
    yield output

  def __call__(*args, **kwargs):
    """ Performs a vasp relaxation. 

        Arguments are the same as `generator`.

        The plane wave basis depends upon cell-shape. Hence, the basis goes out
        of sync with the actual structure during cell-shape relaxation.
        Convergence can only be achieved by restarting the calculation. And
        eventually performing a static calculation.

        It is wastefull to use this functional unless a cell-shape relaxation
        is wanted.

        On top of the usual vasp parameter, this functional accepts the keyword
        agument *relaxation* as a shorthand for *set_relaxation* attribute of
        the vasp functional. If both are passed to this method,
        *set_relaxation* takes precedence.

        If you want to examine the result of each and every vasp calculation,
        use `generator` instead.
    """ 
    for output in self.generator(*args, **kwargs): pass
    return output


# def kpoint_convergence(vasp, structure, outdir="kconv", comm = None, start=1, steps=None, \
#                        offset=(0.5,0.5,0.5), repat=[], tolerance=1e-3, **kwargs):
#   """ Performs a convergence test for kpoints using kpoints.Density object.

#       This is a generator which yields a output extraction object after each
#       vasp calculation:
#       >>> for extract in kpoint_convergence(structure, vasp, outdir)
#       >>>   print extract.total_energy
#       Note that this function will not recompute the result if an appropriate
#       output directory is found. 
#       @note: This method works only if C{vasp.kpoints = integer(??)} is correct.
#       @note: This functor is stateless as long as self and structure can be
#              deepcopied correctly.  

#       @param structure: A structure to relax.
#       @type structure: L{lada.crystal.Structure} 
#       @param vasp: The vasp functional.
#       @type vasp: L{Vasp}
#       @param outdir: Directory where to repatriate files. Default = kconv.
#       @type outdir: str
#       @param comm: Groups of processes for which to call vasp.
#       @type comm: None or boost.mpi.Communicator 
#       @param start: Starting density. Default = 1.
#       @param steps: End density. Default = 1.
#       @param offset: Offset from L{Gamma} of reciprocal mesh. Default = (0.5,0.5,0.5). 
#       @type offset: 3-tuple.
#       @param repat: File to repatriate, other than L{files.minimal}. Default: [].
#       @type repat: list or set
#       @param tolerance: Total energy convergence criteria (per atom). Default: 1e-3. 
#       @type tolerance: float
#   """
#   from copy import deepcopy
#   from math import fabs as abs
#   from os.path import join, exists
#   from os import remove
#   from ..kpoints import Density
#   from .. import files

#   # make this function stateless.
#   vasp = deepcopy(vasp)
#   repat = set(repat).union(files.minimal)
#   tolerance = float(tolerance) * float( len(structure.atoms) )
#   vasp.relaxation.value = "static" # what else could it be?
#   density = deepcopy(start)

#   # keywords arguments cannot include kpoints.
#   if "kpoints" in kwargs:
#     raise SyntaxError,\
#           "Cannot have kpoints as a keyword argument when performing kpoints convergence.\n"

#   # we may want to include restart file to speed up calculations. However,
#   # these files should be deleted afterwards.
#   other_repat = []
#   for file in files.restart: 
#     if file not in repat: other_repat.append(file)
#   repat.union(other_repat)

#   # performs initial calculation.
#   outdirs = ["%s/density:%i" % (outdir, density)]
#   output = vasp( structure, outdirs[-1], comm = comm, repat = repat,\
#                  kpoints=Density(offset, density), **kwargs )
#   # yields output for whatnot
#   yield output 
#   # makes sure we don't accidentally converge to early.
#   oldenergy = -output.total_energy

#   while( abs(oldenergy - output.total_energy) > tolerance ):
#     # restart + new directory
#     vasp.indir = outdirs[-1]
#     density += steps

#     # launch calculations 
#     outdirs.append("%s/density:%i" % (outdir, density))
#     output = vasp( structure, outdirs[-1], comm = comm, repat = repat,\
#                    kpoints=Density(offset, density), **kwargs )
#     # yields output for whatnot.
#     yield output

#     # makes sure we don't accidentally converge to early.
#     oldenergy = output.total_energy

#   # cleanup -- deletes unwanted files from previous output directory
#   for dir in outdirs:
#     for file in other_repat:
#       filename = join(dir, file)
#       if exists(filename): remove(filename)
#   
# def main():
#   import os.path
#   import shutil
#   from boost.mpi import world 
#   from numpy import array as np_array
#   from lada import crystal
#   from lada.vasp import Vasp, Specie

#   structure = crystal.Structure()
#   structure.cell = np_array( [[1.0,0,0],[0,1,0],[0,0,1]] )
#   structure.atoms.append( crystal.StrAtom(np_array([0.0,0,0]), "Rb") )
#   structure.atoms.append( crystal.StrAtom(np_array([0.5,0.5,0.5]), "K") )
#   structure.name = "KRb_nosym"
#   structure.scale = 6

#   K = Specie( "K", "~/AtomicPotentials/pseudos/K_s" )
#   Rb = Specie( "Rb", "~/AtomicPotentials/pseudos/Rb_s" )

#   vasp = Vasp(species=(K, Rb))
#   if os.path.exists(vasp.indir):  shutil.rmtree( vasp.indir )

#   for output in relaxation( structure, vasp, outdir = "test", comm = world ):
#     structure = output.structure
#     print structure.name, structure.energy, "\n", structure

# if __name__ == "__main__":
#   main()

#     

#     
#  
