""" Sub-package containing the functional. """
__docformat__ = "restructuredtext en"
__all__ = ['Functional']
from ..functools import stateless, assign_attributes
from ..misc import add_setter
from extract import Extract
from incar import Incar

class Functional(Incar):
  """ Interface to VASP code. """
  Extract = staticmethod(Extract)
  """ Extraction class. """

  def __init__(self, copyfrom=None, species=None, kpoints=None, **kwargs):
    """ Initializes vasp class. """
    super(Functional, self).__init__(**kwargs)

    self.restart_from_contcar = kwargs.pop('restart_from_contcar', True)
    """ If True and self.CONTCAR exists in directory, will restart from it. """

    # copies values from other functional.
    if copyfrom is not None: 
      self.params.update(copyfrom.params)
      self.special.update(copyfrom.special)
      for key, value in copyfrom.__dict__.iteritems():
        if key in kwargs: continue
        elif key == 'params': continue
        elif key == 'special': continue
        elif hasattr(self, key): setattr(self, key, value)

    self.species = species if species is not None else {}
    """ Species in the system. """
    self.kpoints = kpoints if kpoints is not None \
                   else "0\nM\n4 4 4\n0 0 0"
    """ kpoints for which to perform calculations. """

    if 'program' in kwargs:
      self.program = kwargs['program']
      """ Name/fullpath of vasp program. """

    # sets all known keywords as attributes.
    for key, value in kwargs.iteritems():
      if hasattr(self, key): setattr(self, key, value)

  def __call__( self, structure, outdir=None, comm=None, overwrite=False, 
                mpirun_exe=None, **kwargs):
    """ Calls vasp program. """
    from ..misc import execute_program
    result = None
    for program in self.iter(structure, outdir=outdir, comm=comm, overwrite=overwrite, **kwargs):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False):
        result = program
        continue
      # otherwise, it should yield a Program tuple to execute.
      execute_program(program, append=False, comm=comm, **kwargs)
      result = self.Extract(outdir)
      if not result.success: raise RuntimeError("Vasp failed to execute correctly.")
    return result

  @assign_attributes(ignore=['overwrite'])
  @stateless
  def iter( self, structure, outdir=None, comm=None, overwrite=False, **kwargs ):
    """ Performs a vasp calculation 
     
        If successfull results (see ``extract.Extract.success``) already exist
        in outdir, calculations are not repeated. Instead, an extraction object
        for the stored results are given.

        :param structure:  
            :class:`Structure <lada.crystal.Structure>` structure to compute,
            *unless* a CONTCAR already exists in ``outdir``, in which case this
            parameter is ignored. (This feature can be disabled with the
            keyword/attribute ``restart_from_contcar=False``).
        :param outdir:
            Output directory where the results should be stored.  This
            directory will be checked for restart status, eg whether
            calculations already exist. If None, then results are stored in
            current working directory.
        :param comm:
            Holds arguments for executing VASP externally.
        :param overwrite:
            If True, will overwrite pre-existing results. 
            If False, will check whether a successfull calculation exists. If one does, 
            then does not execute. 
        :param kwargs:
            Any attribute of the VASP instance can be overidden for
            the duration of this call by passing it as keyword argument.  

        :return: Yields an extractor object if a prior successful run exists.
                 Otherwise, yields a tuple object for executing an external
                 program.

        :note: This functor is stateless as long as self and structure can be
               deepcopied correctly.  

        :raise RuntimeError: when computations do not complete.
        :raise IOError: when outdir exists but is not a directory.
    """ 
    from os.path import exists, join
    from numpy import abs
    from numpy.linalg import det
    from ..crystal import specieset
    from ..crystal import read
    from .files import CONTCAR
    from .. import vasp_program
    from ..misc import Program

    # check for pre-existing and successfull run.
    if not overwrite:
      extract = self.Extract(outcar=outdir)
      if extract.success:
        yield extract # in which case, returns extraction object.
        return
    
    # makes functor stateless/reads structure from CONTCAR if requested and appropriate.
    if kwargs.pop("restart_from_contcar", self.restart_from_contcar): 
      path = join(outdir, CONTCAR)
      if exists(path):
        try: contstruct = read.poscar(path, specieset(structure))
        except: pass
        else: structure = contstruct
    if len(structure) == 0: raise ValueError("Structure is empty.")
    if abs(det(structure.cell)) < 1e-8: raise ValueError("Structure with zero volume.")
    if abs(structure.scale) < 1e-8: raise ValueError("Structure with null scale.")

    self.pullup(structure, outdir, comm)
    program = getattr(self, 'program', vasp_program)
    yield Program(program, [], outdir, 'stdout', 'stderr')
    self.bringdown(outdir, structure)

  def pullup(self, structure, outdir, comm):
    """ Creates all input files necessary to run results.

        Performs the following actions.

        - Writes POSCAR file.
        - Writes INCAR file.
        - Writes KPOINTS file.
        - Creates POTCAR file
        - Saves pickle of self.
    """
    import cPickle
    from . import files
    from ..crystal import specieset, write
    from ..misc.changedir import Changedir

    with Changedir(outdir) as tmpdir:
      # creates poscar file. Might be overwriten by restart.
      write.poscar(structure)

      # creates incar file. Changedir makes sure that any calculations done to
      # obtain incar will happen in the right directory.
      incar_lines = self.incar_lines(structure=structure, vasp=self, comm=comm)

      # creates INCAR file. Note that POSCAR file might be overwritten here by Restart.
      with open(files.INCAR, "w") as incar_file: 
        incar_file.writelines(incar_lines)
  
      # creates kpoints file
      with open(files.KPOINTS, "w") as kp_file: 
        self.write_kpoints(kp_file)
  
      # creates POTCAR file
      with open(files.POTCAR, 'w') as potcar:
        for s in specieset(structure):
          potcar.writelines( self.species[s].read_potcar() )
    
      with open(files.FUNCCAR, 'w') as file:
        cPickle.dump(self, file)
    

  def bringdown(self, directory, structure):
     """ Copies contcar to outcar. """
     from . import files
     from ..misc.changedir import Changedir

     # Appends INCAR and CONTCAR to OUTCAR:
     with Changedir(directory) as pwd:
       with open(files.OUTCAR, 'a') as outcar:
         outcar.write('\n################ CONTCAR ################\n')
         with open(files.CONTCAR, 'r') as contcar: outcar.write(contcar.read())
         outcar.write('\n################ END CONTCAR ################\n')
         outcar.write('\n################ INCAR ################\n')
         with open(files.INCAR, 'r') as incar: outcar.write(incar.read())
         outcar.write('\n################ END INCAR ################\n')
         outcar.write('\n################ INITIAL STRUCTURE ################\n')
         outcar.write(repr(structure))
         outcar.write('\n################ END INITIAL STRUCTURE ################\n')

    
  def write_kpoints(self, file):
    """ Writes kpoints to a stream. """
    if isinstance(self.kpoints, str): file.write(self.kpoints)
    elif hasattr(self.kpoints, "__call__"): file.write(self.kpoints(self))
    else: # numpy array or such.
      file.write("Explicit list of kpoints.\n{0}\nCartesian\n".format(len(self.kpoints)))
      for kpoint in self.kpoints:
        file.write("{0[0]} {0[1]} {0[2]} {1}\n".format(kpoint, 1 if len(kpoint) == 3 else kpoint[3]))

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
      if value is None: continue
      # if a special parameter, then is non-default.
      if name in params: string += "functional.%s = %s\n" % (name, repr(value))
      else:
        string += "functional.add_param = \"%s\", %s\n" % (name, repr(value))
        module = value.__class__.__module__ 
        classname = value.__class__.__name__ 
        if module in modules: modules[module].append(classname)
        else: modules[module] = [classname]
    for name, value in self.special.items():
      if value.value is None: continue
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
    if not self.restart_from_contcar: 
      string += "functional.restart_from_contcar = False\n"

    # adds user modules above repr string.
    header = ""
    for name in sorted(modules.keys()):
      mods = list(set(modules[name]))
      header += "from %s import %s" % (name, mods[0])
      for v in mods[1:]: header += ", %s" % (v)
      header += "\n"
    return header + string

  @add_setter
  def add_specie(self, args):
    """ Adds a specie to current functional. 
     
        The argument is a tuple containing the following.

        - Symbol (str).
        - Directory where POTCAR resides (str).
        - List of U parameters (optional, see module vasp.specie).
        - Maximum (or minimum) oxidation state (optional, int).
        - ... Any other argument in order of `vasp.specie.Specie.__init__`.
    """
    from .specie import Specie
    assert len(args) > 1, ValueError("Too few arguments.")
    self.species[args[0]] = Specie(*args[1:])
    
  def __setstate__(self, args):
    """ Sets state from pickle.

        Takes care of older pickle versions.
    """
    super(Functional, self).__setstate__(args)
    for key, value in self.__class__().__dict__.iteritems():
       if not hasattr(self, key): setattr(self, key, value)

