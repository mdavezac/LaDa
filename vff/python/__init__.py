""" Valence Force Field functional for zinc-blende materials """

# imports C++ extension
from ..opt.decorators import add_setter, broadcast_result, make_cached
from ..opt import AbstractExtractBase

def _get_script_text(file, name):
  string = "# " + name + " definition."
  for line in file:
    if line.find(string) != -1: break;
  lines = ""
  string = "# End of " + name.lower() + " definition."
  for line in file:
    if line.find(string) != -1: break;
    lines += line
  return lines

def exec_input(filepath = "input.py", namespace = None):
  """ Executes an input script including namespace for escan/vff. """ 
  from ..opt import exec_input

  dictionary = { "Vff": Vff }
  if namespace != None: dictionary.update(namespace)
  return exec_input(filepath, dictionary)

def read_input(filepath = "input.py", namespace = None, name=None):
  """ Reads an input file including namespace for escan/vff. """ 
  from ..opt import read_input

  dictionary = { "Vff": Vff }
  if namespace != None:
    dictionary.update(namespace)
    if name == None and hasattr(namespace, '__name__'): namee = namespace.__name__
  return read_input(filepath, dictionary)

class Extract(AbstractExtractBase):
  """ Extracts vff results from output file. """
  def __init__(self, directory = None, comm = None, vff = None):
    """ Initializes the extraction class. """
    super(Extract, self).__init__(directory=directory, comm=comm)

    self.OUTCAR = vff.OUTCAR if vff != None else Vff().OUTCAR
    """ Filename of the OUTCAR file from VASP.
     
        Data will be read from directory/OUTCAR. 
    """
    self.FUNCCAR = vff._FUNCCAR if vff != None else Vff()._FUNCCAR
    """ Pickle filename for the functional. """

  def __outcar__(self):
    """ Path to OUTCAR. """
    from os.path import exists, join, isfile
    path = join(self.directory, self.OUTCAR)
    assert exists(path), IOError('Path {0} does not exist.'.format(path))
    assert isfile(path), IOError('Path {0} is not a file.'.format(path))
    return open(path, 'r')

    
  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for VFF success.
        
        At this point, checks for files and 
    """
    try:
      with self.__outcar__() as file:
        for line in file:
          if line.find("# Computed VFF in:") != -1: return True
    except: pass
    return False
 
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def structure(self):
    """ Greps structure from self.L{OUTCAR}. """
    with self.__outcar__() as file:
      # find start of calculations.
      for line in file:
        if line.find("# Result of VFF calculations.") != -1: break
      script = _get_script_text(file, "Structure")
    return exec_input(script).structure

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def input_structure(self):
    """ Greps input structure from self.L{OUTCAR}. """
    with self.__outcar__() as file:
      # find start of calculations.
      for line in file:
        if line.find("# Input Structure.") != -1: break
      script = _get_script_text(file, "Structure")
    return exec_input(script).structure

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def energy(self):
    """ Greps energy from self.L{OUTCAR}. """
    with self.__outcar__() as file:
      # find start of calculations.
      for line in file:
        if line.find("# Result of VFF calculations.") != -1: break
      for line in file:
        if line.find("# Structure definition.") != -1: break;
      for line in file:
        if line.find("structure.energy") != -1:
          return float(line.split()[2])

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def stress(self):
    """ Greps stress from self.L{OUTCAR}. """
    result = []
    with self.__outcar__() as file:
      # find start of calculations.
      for line in file:
        if line.find("# Result of VFF calculations.") != -1: break
      lines = []
      for line in file:
        if line.find("stress: ") != -1: lines.append(line); break
      for i, line in zip(range(2), file):
        lines.append(line)
      for line in lines:
        line = line.replace("stress:",'')
        line = line.replace('(','')
        line = line.replace(')','')
        line = line.replace(',','')
        line = line.replace('\\','')
        result.append( [ float(u) for u in line.split() ] )
    return array(result, dtype="float64")

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def lattice(self):
    """ Greps lattice from self.L{OUTCAR}. """
    from lada.crystal import Lattice
    with self.__outcar__() as file: script = _get_script_text(file, "Lattice")
    return exec_input(script).lattice

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def minimizer(self):
    """ Greps minimizer from self.L{OUTCAR}. """
    from ..minimizer import Minimizer
    with self.__outcar__() as file: script = _get_script_text(file, "Minimizer")
    return exec_input(script).minimizer

  @property
  @make_cached
  def vff(self):
    """ Greps vff functional from self.L{OUTCAR}. """
    from os.path import exists, join
    from cPickle import load
    from numpy import array
    from . import Vff

    # tries to read from pickle.
    path = self.FUNCCAR
    if len(self.directory): path = join(self.directory, self.FUNCCAR)
    if exists(path):
      try:
        with open(path, "r") as file: result = load(file)
      except: pass 
      else: return result

    # tries to read from outcar.
    @broadcast_result(attr=True, which=0) # easier to broadcast this way.
    def get_vff(this):
      with self.__outcar__() as file: return _get_script_text(file, "Vff")

    local_dict = {"lattice": self.lattice, "array": array, "minimizer": self.minimizer, "Vff": Vff}
    input = exec_input(get_vff(self))
    return input,vff_functional if "vff_functional" in input else input.functional

  def write_escan_input(self, filepath, structure = None):
    """ Prints escan input to file. """
    from os.path import expanduser
    if structure == None: structure = self.structure

    # Saves global lattice if set.
    old_lattice = None
    try: old_lattice = structure.lattice
    except RuntimeError: pass
    self.lattice.set_as_crystal_lattice()

    functional = self.vff._create_functional()
    functional.print_escan_input(expanduser(filepath), structure)

    if old_lattice != None: old_lattice.set_as_crystal_lattice()


class Vff(object): 
  """ Valence Force Field functional for zinc-blende materials """
  def __init__(self, workdir = None):
    """ Initializes a valence force field functional. """
    from ..crystal.binary import zinc_blende
    from ..minimizer import Minimizer

    super(Vff, self).__init__()
    self.lattice = zinc_blende()
    """ Lattice for which to perform calculations.
    
        In practice, zinc-blende lattices can be defined for any number of parameterizations.
        The one chosen here is such that atomic site 0 is at the origin of the
        unit cell, and atomic site 1 at 0.25, 0.25, 0.25. The scale is
        arbitrary at this point, and set to 4.5. The unit cell is [[0,0.5,0.5],
        [0.5,0,0.5], [0.5,0.5,0]]. 
    """

    self.relax = True
    """ Whether to relax strain or just get energy. """
    self.direction = None
    """ Epitaxial direction if any.
    
        If set, it should be a vector. In that case, relaxation will only
        happen in that direction. Otherwise, everything will be relaxed. 
    """
    self.minimizer = Minimizer()
    """ Minimizer to use with VFF """
    self.bonds = {}
    """ Bond parameters. """
    self.angles = {}
    """ Angle parameters. """
    self.OUTCAR = "vff_out" 
    """ File where to redirect output. """
    self.ERRCAR = "vff_err" 
    """ File where to redirect errors. """
    from cPickle import dump
    self._FUNCCAR = "VFFCAR" 
    """ Pickle file to which functional is saved. """
    self._workdir = workdir
    """ Private reference to the working directory. """
    self.print_from_all = False
    """ If True, each node will print. """

  def _get_workdir(self): return self._workdir
  def _set_workdir(self, workdir):
    from os.path import abspath, expanduser
    self._workdir = abspath(expanduser(workdir)) if workdir != None else None
  workdir = property( _get_workdir, _set_workdir, 
                      """ Working directory where calculations are performed. 
                      
                          Absolute path is determined when workdir is set.
                      """ )

  @add_setter
  def add_bond(self, args):
    """ Adds/Modifies the bond parameter dictionary.
    
        >>> vff.add_bond = "A", "B", [d0, a2, a3, a4, a5, a6]
        @param args: - the first element of arg is one endpoint of the bond (say "In"), 
                     - the second element is another endpoint (say "N"),
                     - the last element is a sequence with at most 5 elements
                       defining the bond-stretching parameters (order 2 through 6).
        @type args: "type1", "type2", sequence
    """
    assert len(args) == 3, RuntimeError("Bonds should be set with 3 parameters: %s." % (args))
    assert args[0] in [ u for site in self.lattice.sites for u in site.type ],\
           RuntimeError( "Type %s is not in lattice." % (args[0]))
    assert args[1] in [ u for site in self.lattice.sites for u in site.type ],\
           RuntimeError( "Type %s is not in lattice." % (args[1]))

    name = "%s-%s" % (args[0], args[1])
    if args[0] > args[1]: name = "%s-%s" % (args[1], args[0])

    params = [u for u in args[2]]
    assert len(params) < 7 and len(params) > 0, RuntimeError("To few or too many parameters: %s." % (params))
    if name in self.bonds: # replaces none value with known value
      for old, new in zip(self.bonds[name], params):
        if old == None: old = new
    if len(params) < 6: params.extend( 0 for u in range(6-len(params)) )
    
    self.bonds[name] = params

  def set_bond(self, A, B, d0 = None, a2 = None, a3 = None, a4 = None, a5 = None, a6 = None):
    """ Adds/Modifies the bond parameter dictionary.
    
        @param A: Endpoint specie.
        @type A: str
        @param B: Other endpoint species. 
        @type B: str
        @param d0: ideal bond-length.
        @param a2: 2nd order bond-stretching parameter.
        @param a3: 3rd order bond-stretching parameter.
        @param a4: 4rd order bond-stretching parameter.
        @param a5: 5rd order bond-stretching parameter.
        @param a6: 6rd order bond-stretching parameter.
    """
    self.add_bond = A, B, [d0, a2, a3, a4, a5, a6]




  @add_setter
  def add_angle(self, args):
    """ Adds/Modifies the angle parameter dictionary.
    
        >>> vff.add_angle = "A", "B", "C", [gamma, sigma, a2, a3, a4, a5, a6]
        @param args: - the first element of arg is one endpoint of the angle (say "In"), 
                     - the second element is middle of the angle (say "N"),
                     - the third element is the other endpoint (say "Ga"),
                     - the last element is a sequence with at most 7 elements
                       defining the bond-angle parameters. The first is the
                       gamma parameter, the second the sigma, and the others
                       the bond bending proper.
        @type args: "type1", "type2", "type3", sequence
    """
    assert len(args) == 4, RuntimeError("Angle should be set with 4 parameters: %s." % (args))
    assert args[0] in [ u for site in self.lattice.sites for u in site.type ],\
           RuntimeError( "Type %s is not in lattice." % (args[0]))
    assert args[1] in [ u for site in self.lattice.sites for u in site.type ],\
           RuntimeError( "Type %s is not in lattice." % (args[1]))
    assert args[2] in [ u for site in self.lattice.sites for u in site.type ],\
           RuntimeError( "Type %s is not in lattice." % (args[2]))

    name = "%s-%s-%s" % (args[0], args[1], args[2])
    if args[0] > args[2]: name = "%s-%s-%s" % (args[2], args[1], args[0])

    params = [u for u in args[3]]
    assert len(params) <= 7 and len(params) > 0,\
           RuntimeError("To few or too many parameters: %s." % (params))
    if name in self.angles: # replaces none value with known value
      for old, new in zip(self.angles[name], params):
        if old == None: old = new
    if str(params[0]).lower()[:3] == "tet": params[0] = -1e0/3e0
    if len(params) < 7: params.extend( 0e0 for u in range(7-len(params)) )

    self.angles[name] = params
    
  def set_angle(self, A, B, C, gamma = None, sigma = None, a2 = None,\
                      a3 = None, a4 = None, a5 = None, a6 = None):
    """ Adds/Modifies the angle parameter dictionary.
    
        @param A: Endpoint specie.
        @type A: str
        @param B: Middle specie.
        @type B: str
        @param C: Endpoint specie.
        @type C: str
        @param gamma: S{gamma} parameter. Can be a number or a string starting with "tet".
        @param sigma: S{sigma} parameter.
        @param a2: 2nd order bond-angle parameter.
        @param a3: 3rd order bond-angle parameter.
        @param a4: 4rd order bond-angle parameter.
        @param a5: 5rd order bond-angle parameter.
        @param a6: 6rd order bond-angle parameter.
    """
    self.add_angle = A,B, C, [gamma, sigma, a2, a3, a4, a5, a6]


  def __repr__(self):
    result  = repr(self.lattice)
    result += repr(self.minimizer)
    result += "# Vff definition.\n"
    result += "functional = %s()\n" % (self.__class__.__name__)
    result += "functional.minimizer = minimizer\n"
    result += "functional.lattice = lattice\n"
    result += "functional.print_from_all = {0}\n".format(repr(self.print_from_all))
    result += "functional.OUTCAR = '%s'\n" % (self.OUTCAR)
    result += "functional.ERRCAR = '%s'\n" % (self.ERRCAR)
    if self.direction == None or hasattr(self.direction, "__len__"):
      result += "functional.direction = %s\n" % (repr(self.direction))
    else: result += "functional.direction = %i\n" % (int(self.direction))
    result += "functional.relax = %s\n" % ("True" if self.relax else "False")
    for name, params in self.bonds.items():
      A, B = name.split("-")
      result += "functional.add_bond = '%s', '%s', (%e" % (A, B, params[0])
      for i, u in enumerate(params[1:]): 
        if sum([abs(j) for j in params[i+1:]]) / float(len(params)-i-1) < 1e-12: break;
        result += ", %e" % (u)
      result += ")\n"
    for name, params in self.angles.items():
      A, B, C = name.split("-")
      result += "functional.add_angle = '%s', '%s', '%s', " % (A, B, C)
      if abs(params[0] + 1e0/3e0) < 1e-12: result += "('tet'"
      else: result += "(%e" % (params[0])
      for i, u in enumerate(params[1:]):
        if sum([abs(j) for j in params[i+1:]]) / float(len(params)-i-1) < 1e-12: break;
        result += ", %e" % (u)
      result += ")\n"
    result += "# End of vff definition.\n"


    module = self.__class__.__module__ 
    classname = self.__class__.__name__ 
    header = "from %s import %s\n" % (module, classname)
    module = self.minimizer.__class__.__module__ 
    classname = self.minimizer.__class__.__name__ 
    header += "from %s import %s\n" % (module, classname)
    return header + result

  def _cout(self, comm):
    """ Creates output name. """
    if self.OUTCAR == None: return "/dev/null"
    if comm == None:   return self.OUTCAR
    if comm.rank == 0: return self.OUTCAR
    return self.OUTCAR + "." + str(comm.rank) if self.print_from_all else "/dev/null"
  def _cerr(self, comm):
    """ Creates error name. """
    if self.ERRCAR == None: return "/dev/null"
    if comm == None:   return self.ERRCAR
    if comm.rank == 0: return self.ERRCAR
    return self.ERRCAR + "." + str(comm.rank) if self.print_from_all else "/dev/null"

  def __call__(self, structure, outdir = None, comm = None, overwrite=False, **kwargs):
    """ Performs calculation """
    import time
    from cPickle import dump
    from copy import deepcopy
    from os import getcwd
    from os.path import exists, isdir, abspath, expanduser
    from .. import lada_with_mpi
    from ..opt.changedir import Changedir

    if lada_with_mpi and comm == None: 
      from boost.mpi import world
      comm = world
    is_mpi  = False if comm == None else comm.size > 1
    is_root = comm.rank == 0 if is_mpi else True

    # bull shit. 
    assert len(self.lattice.sites) == 2, RuntimeError("Lattice is not zinc-blend")
    assert len(self.lattice.sites[0].type) > 0,\
           RuntimeError("No atomic species given in lattice site 0.")
    assert len(self.lattice.sites[1].type) > 0,\
           RuntimeError("No atomic species given in lattice site 1.")
    assert len(self.lattice.sites[1].type) < 3,\
           RuntimeError("Lattice is not binary, pseudo-binary or quaternary.")
    assert len(self.lattice.sites[0].type) < 3,\
           RuntimeError("Lattice is not binary, pseudo-binary or quaternary.")
    assert len(self.lattice.sites[1].type) < 3,\
           RuntimeError("Lattice is not binary, pseudo-binary or quaternary.")
    assert len(self.lattice.sites[0].type) != 1 or  len(self.lattice.sites[1].type) != 2

    timing = time.time() 
    local_time = time.localtime() 

    # gets absolute path.
    outdir = abspath(expanduser(outdir)) if outdir != None else getcwd()

    # make this functor stateless.
    this      = deepcopy(self)
    assert abs(structure.scale) > 1e-12 or abs(self.lattice.scale) > 1e-12,\
           RuntimeError("Scales in input structure and lattice are both zero.")
    if abs(structure.scale) < 1e-12: structure.scale = self.lattice.scale

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to use for calculations launch. 
    # For added simplicity, minimizer attributes are also directly available.
    for key in kwargs.keys():
      if hasattr(this, key): setattr(this, key, kwargs[key])
      elif hasattr(this.minimizer, key): setattr(this.minimizer, key, kwargs[key])
      else: raise NameError( "%s attribute unknown of vff or vff.minimizer." % (key) )


    if comm == None: comm = world

    # checks if outdir contains a (wanted) successful run.
    does_exist = exists(outdir) if is_root else None
    if is_mpi: 
      from boost.mpi import broadcast
      does_exist, overwrite = broadcast(comm, (does_exist, overwrite), 0)
    if does_exist and not overwrite:
      extract = Extract(comm = comm, directory = outdir, vff = this)
      if extract.success: return extract # in which case, returns extraction object.
    if is_mpi: comm.barrier() # makes sure directory is not created by other proc!
    
    with Changedir(outdir, comm = comm) as current_dir:
      # redirects C/C++/fortran streams
      cout, cerr = this._cout(comm), this._cerr(comm)
      # write stuff to output file
      with open(cout, "w") as file:
        # now for some input variables
        print >> file, "# VFF calculation on ", time.strftime("%m/%d/%y", local_time),\
                       " at ", time.strftime("%I:%M:%S %p", local_time)
        print >> file, "# Input Structure.\n{0}\n".format(repr(structure))
        if len(structure.name) != 0: print "# Structure named ", structure.name 
        print >> file, repr(this)
        print >> file, "# Performing VFF calculations. "
        # then calculations

      # Saves FUNCCAR.
      if is_root:
        with open(this._FUNCCAR, "w") as file: dump(this, file)

      result, stress = this._run(structure, comm)
      # must close/reopen redirection context, otherwise it seem that they 
      # are closed anyhow on exiting from the Cpp function call. The following context
      with open(cout, "r") as file: lines = file.readlines()
      with open(cout, "w") as file:
        file.writelines(lines)
        # finally, the dam results.
        file.write("# Result of VFF calculations.\n")
        file.write(repr(result))
        if stress != None: 
          file.write( "\nstress: ({0[0]}, {0[1]}, {0[2]}),\\"\
                      "\n        ({1[0]}, {1[1]}, {1[2]}),\\"\
                      "\n        ({2[0]}, {2[1]}, {2[2]})\n"\
                      .format(stress[0,:], stress[1,:], stress[2,:]))
        else: file.write("\nstress not computed.\n")
        timing = time.time() - timing
        hour = int(float(timing/3600e0))
        minute = int(float((timing - hour*3600)/60e0))
        second = (timing - hour*3600-minute*60)
        file.write("# Computed VFF in: %i:%i:%f.\n"  % (hour, minute, second))


    # checks if result was successful
    extract = Extract(comm = comm, directory = outdir, vff = this)
    assert extract.success, RuntimeError("VFF calculation did not complete in %s.\n" % (outdir))

    return extract


  def _create_functional(self):
    """ Creates the vff functional using cpp extension. """
    from tempfile import NamedTemporaryFile
    from os import remove
    from _vff import Vff, LayeredVff
    from sys import exit

    if hasattr(self.direction, "__len__"):
      functional = Layered()
      functional.direction = self.direction
    else: functional = Vff()
    # minimizer variants are somewhat difficult to expose...
    self.minimizer._copy_to_cpp(functional._minimizer)
    # ... done jumping through hoops.
    for name, params in self.bonds.items():
      bond = functional._get_bond(name.split('-'))
      bond.length = params[0]
      bond.length, bond.alphas[:] = params[0], params[1:]
    for name, params in self.angles.items():
      angle = functional._get_angle(name.split('-'))
      angle.gamma, angle.sigma, angle.betas[:] = params[0], params[1], params[2:]

    return functional
     
  def _run(self, structure, comm):
    """ Performs actual calculation. """
    from copy import deepcopy
    from _vff import Vff, LayeredVff
    from .. import lada_with_mpi
    from ..opt import redirect_all

    # Saves global lattice if set.
    old_lattice = None
    try: old_lattice = structure.lattice
    except RuntimeError: pass
    self.lattice.set_as_crystal_lattice()
    # if direction is not None and not an array, then should be a combination
    # of FreezeCell integers.
    if self.direction != None and not hasattr(self.direction, "__len__"):
      oldfreeze = structure.freeze
      structure.freeze = self.direction
    
    cout, cerr = self._cout(comm), self._cerr(comm)
    with open(cerr, "w") as file: pass # file has not yet been opened
    with redirect_all(output=cout, error=cerr, append="True") as oestream:
      functional = self._create_functional()
    # now performs call
    functional.init(structure, dotree=True)
    functional.check_input()
    
    with redirect_all(output=cout, error=cerr, append="True") as oestream:
      result, stress = functional(comm, relax=self.relax) if lada_with_mpi\
                       else functional(relax=self.relax)
    
    # unsets lattice.
    if old_lattice != None: old_lattice.set_as_crystal_lattice()
    if self.direction != None and not hasattr(self.direction, "__len__"):
      structure.freeze = oldfreeze
    return result, stress
