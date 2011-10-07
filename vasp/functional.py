""" Sub-package containing the functional. """
__docformat__ = "restructuredtext en"
__all__ = ['Functional']
from launch import Launch
from extract import Extract

class Functional(Launch):
  """ Interface to VASP code.
     
      The strength of this interface is that combines vasp calculations, result
      caching, and data retrieval together. 
      
      A vasp run is parameterized using Incar class defined in incar.py.
      It is launched using the Launch class from launch.py class. More
      specifically, VASP derives from Launch wich derives from Incar.  The
      results of a successful run is cached in the self.outdir directory.
      After being launched an object is returned which can extract output data
      from the files in this directory.

      One of the strength of this class is that since results are cached in the
      outdir directory, successful calculations are never runned twice.
      This allows us to use the same scripts for generating and retrieving
      data. 
  """
  Extract = staticmethod(Extract)
  """ Extraction class. """

  def __init__(self, vasp=None, **kwargs):
    """ Initializes vasp class. """
    super(Functional, self).__init__(self, **kwargs)

    self.restart_from_contcar = kwargs.pop('restart_from_contcar', True)
    """ If True and self.CONTCAR exists in directory, will restart from it. """

    # copies values from other functional.
    if vasp != None: 
      self.params.update(vasp.params)
      self.special.update(vasp.special)
      for key, value in vasp.__dict__.iteritems():
        if key in kwargs: continue
        elif key == 'params': continue
        elif key == 'special': continue
        elif hasattr(self, key): setattr(self, key, value)

    # sets all known keywords as attributes.
    for key, value in kwargs.iteritems():
      if hasattr(self, key): setattr(self, key, value)

  def __call__( self, structure, outdir = None, comm = None, repat = None,\
                overwrite=False, keep_tempdir=False, minversion=0, **kwargs ):
    """ Performs a vasp calculation 
     
        :Parameters:
	  structure :  `lada.crystal.Structure`
            the structure to compute, *unless* a CONTCAR already exists in
            ``outdir``, in which case this parameter is ignored. (This feature
            can be disabled with the keyword/attribute
            ``restart_from_contcar=False``).
	  outdir : str or None
            Output directory where the results should be stored.  This
            directory will be checked for restart status, eg whether
            calculations already exist. If None, then results are stored in
            current working directory.
          comm : `mpi.Communicator`
            Calculations are performed over this MPI communicator.  
          repat : list 
            list of files to save, above and beyond `files.minimal`. This is
            only used when `Functional.inplace` is False.
	  overwrite : boolean 
            Whether to overwrite pre-existing results, eg does not check for
            restart status.
	  kwargs 
            Any attribute of the VASP instance can be overidden for
            the duration of this call by passing it as keyword argument.  

        :return: An `extract.Extract` object from which results can be obtained.

        If successfull results (see ``extract.Extract.success``) already exist
        in outdir, calculations are not repeated. Instead, an extraction object
        for the stored results are given.

        :note: This functor is stateless as long as self and structure can be
               deepcopied correctly.  

        :raise RuntimeError: when computations do not complete.
        :raise IOError: when outdir exists but is not a directory.
    """ 
    from copy import deepcopy
    from os import getcwd
    from os.path import exists, join
    from numpy import abs
    from numpy.linalg import det
    from ..crystal import specie_list, read_poscar
    from ..opt import RelativeDirectory
    from ..mpi import Communicator
    from .files import CONTCAR

    # make this functor stateless.
    this      = deepcopy(self)
    outdir    = getcwd() if outdir == None else RelativeDirectory(outdir).path
    repat     = deepcopy(repat)  if repat != None else []
    norun     = kwargs.pop("norun", False)
    # makes functor stateless/reads structure from CONTCAR if requested and appropriate.
    if "external" in kwargs: this.launch_as_library = not kwargs.pop("external")
    if kwargs.pop("restart_from_contcar", self.restart_from_contcar): 
      path = join(outdir, CONTCAR)
      if exists(path):
        try:
          old_structure = structure
          structure = read_poscar(specie_list(structure), path)
          assert len(structure.atoms) > 0
          assert abs(det(structure.cell)) > 1e-8
        except: structure = deepcopy(old_structure)
    else: structure = deepcopy(structure) 

    comm = Communicator(comm)

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to be changed before launch. 
    for key, value in kwargs.items():
      # direct attributes.
      if hasattr(this, key): setattr(this, key, value)
      # properties attributes.
      elif hasattr(this.__class__, key): setattr(this, key, value)
      else: raise ValueError("Unkwown keyword argument to vasp: %s=%s" % (key, value))

    # Checks for previous run, or deletes previous run if requested.
    if not overwrite:
      extract = self.Extract(comm = comm, outcar = outdir)
      if extract.success: return extract # in which case, returns extraction object.
    comm.barrier() # sync all procs.
    
    # Otherwise, performs calculation by calling base class functor.
    super(Functional, this).__call__( structure=structure, outdir=outdir,\
                                      repat=repat, comm=comm, norun=norun, \
                                      keep_tempdir=keep_tempdir, minversion=minversion )
    
    # checks if result was successful
    extract = self.Extract(comm = comm, outcar = outdir)
    if not norun:
      assert extract.success, RuntimeError("VASP calculation did not complete in %s.\n" % (outdir))

    return extract

  def __repr__(self):
    """ Returns a python script representing this object. """
    from .incar._params import SpecialVaspParam
    string = "functional = %s()\n" % (self.__class__.__name__)

    # creates a default vasp instance to compare to.
    compare = self.__class__()
    params = compare.params.keys()

    # will hold classes from modules.
    modules = {}
    modules[self.__class__.__module__] = [self.__class__.__name__]
    # now go through vasp parameters and print them out.
    for name, value in self.params.items():
      if value == None: continue
      # if a special parameter, then is non-default.
      if name in params: string += "functional.%s = %s\n" % (name, repr(value))
      else:
        string += "functional.add_param = \"%s\", %s\n" % (name, repr(value))
        module = value.__class__.__module__ 
        classname = value.__class__.__name__ 
        if module in modules: modules[module].append(classname)
        else: modules[module] = [classname]
    for name, value in self.special.items():
      if value.value == None: continue
      assert isinstance(value, SpecialVaspParam)
      string += "functional.{0} = {1}\n".format(name, repr(value))
      module = value.__class__.__module__ 
      classname = value.__class__.__name__ 
      if module in modules: modules[module].append(classname)
      else: modules[module] = [classname]
    # adds kpoints
    string += "functional.kpoints = %s\n" % (repr(self.kpoints))
    if hasattr(self.kpoints, "__call__"):
      # checks for user module.
      module = self.kpoints.__class__.__module__ 
      classname = self.kpoints.__class__.__name__ 
      if module in modules: modules[module].append(classname)
      else: modules[module] = [classname]
    # adds species.
    for name, specie in self.species.items():
      string += "functional.species['{0}'] = {1}\n".format(name, repr(specie))
      module = specie.__class__.__module__ 
      classname = specie.__class__.__name__ 
      if module in modules: modules[module].append(classname)
      else: modules[module] = [classname]
    if not self.inplace: 
      string += "functional.inplace = False\n"
      string += "functional.workdir = \"%s\"\n" % (self._workdir.repr())
    if not self.restart_from_contcar: 
      string += "functional.restart_from_contcar = False\n"
    if self.print_from_all: 
      string += "functional.print_from_all = True\n"

    # adds user modules above repr string.
    header = ""
    for name in sorted(modules.keys()):
      mods = list(set(modules[name]))
      header += "from %s import %s" % (name, mods[0])
      for v in mods[1:]: header += ", %s" % (v)
      header += "\n"
    return header + string
