""" Valence Force Field functional for zinc-blende materials """

# imports C++ extension
from ..opt.decorators import add_setter, broadcast_result, make_cached

def _is_in_sync(comm, which = [0]):
  from boost.mpi import broadcast
  if comm == None: return 
  print "sync ", comm.rank, which[0]
  which[0] += 1
  m = broadcast(comm, "666" if comm.rank == 0 else None, 0)
  return m == "666"

def zinc_blend_lattice():
  """ Defines a default zinc-blende lattice (InGaNP). """
  from ..crystal import Lattice, Site
  from numpy import array

  lattice = Lattice()
  lattice.set_cell = (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0)
  lattice.scale = 1e0
  lattice.add_site = (0,0,0)
  lattice.add_site = (0.25,0.25,0.25)
  return lattice

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

class Extract(object):
  def __init__(self, directory = None, comm = None, vff = None):
    super(Extract, self).__init__()
    self.directory = directory if directory != None else getcwd()
    """ Directory where to check for output. """
    self.comm = comm
    """ Mpi group communicator. """
    self.OUTCAR = vff.OUTCAR if vff != None else Vff().OUTCAR
    """ Filename of the OUTCAR file from VASP.
     
        Data will be read from directory/OUTCAR. 
    """

    
  def _get_directory(self):
    """ Directory with VASP output files """
    return self._directory
  def _set_directory(self, dir):
    from os.path import abspath, expanduser
    dir = abspath(expanduser(dir))
    if hasattr(self, "_directory"): 
      if dir != self._directory: self.uncache()
    self._directory = dir
  directory = property(_get_directory, _set_directory)


  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for VFF success.
        
        At this point, checks for files and 
    """
    from os.path import exists, join
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    if not exists(path): return False

    with open(path, "r") as file:
      for line in file:
        if line.find("# Computed VFF in:") != -1: return True
    return False
 
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def structure(self):
    """ Greps structure from self.L{OUTCAR}. """
    from os.path import exists, join
    from ..crystal import Structure
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))

    with open(path, "r") as file:
      # find start of calculations.
      for line in file:
        if line.find("# Result of VFF calculations.") != -1: break
      script = _get_script_text(file, "Structure")
    local_dict = {"Structure": Structure}
    exec script in globals(), local_dict
    return local_dict["structure"]


  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def energy(self):
    """ Greps energy from self.L{OUTCAR}. """
    from os.path import exists, join
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))

    with open(path, "r") as file:
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
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))

    result = []
    with open(path, "r") as file:
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
    from os.path import exists, join
    from lada.crystal import Lattice

    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file: script = _get_script_text(file, "Lattice")
    local_dict = {"Lattice": Lattice}
    exec script in globals(), local_dict
    return local_dict["lattice"]

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def minimizer(self):
    """ Greps minimizer from self.L{OUTCAR}. """
    from os.path import exists, join
    from ..minimizer import Minimizer
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file: script = _get_script_text(file, "Minimizer")
    local_dict = {"Minimizer": Minimizer}
    exec script in globals(), local_dict
    return local_dict["minimizer"]

  @property
  @make_cached
  def vff(self):
    """ Greps vff functional from self.L{OUTCAR}. """
    from os.path import exists, join
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))

    @broadcast_result(attr=True, which=0) # easier to broadcast this way.
    def get_vff(this):
      with open(path, "r") as file: return _get_script_text(file, "Vff")

    local_dict = {"lattice": self.lattice, "minimizer": self.minimizer, "Vff": Vff}
    exec get_vff(self) in globals(), local_dict
    return local_dict["vff"]

  def write_escan_input(self, filepath, structure = None):
    """ Prints escan input to file. """
    if structure == None: structure = self.structure

    # Saves global lattice if set.
    old_lattice = None
    try: old_lattice = structure.lattice
    except RuntimeError: pass
    self.lattice.set_as_crystal_lattice()

    functional = self.vff._create_functional(structure, self.comm)
    functional.print_escan_input(filepath, structure)

    if old_lattice != None: old_lattice.set_as_crystal_lattice()

  def solo(self):
    """ Extraction on a single process.

        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. C{self.L{solo}()} returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    from copy import deepcopy
    
    if self.comm == None: return self
    comm = self.comm 
    self.comm = None
    copy = deepcopy(self)
    self.comm = comm
    return copy


class Vff(object): 
  """ Valence Force Field functional for zinc-blende materials """
  def __init__(self, workdir = None):
    """ Initializes a valence force field functional. """
    from ..minimizer import Minimizer

    super(Vff, self).__init__()
    self.lattice = zinc_blend_lattice()
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
    self._workdir = workdir
    """ Private reference to the working directory. """

  def _get_workdir(self): return self._workdir
  def _set_workdir(self, workdir):
    from os.path import abspath, expanduser
    self._workdir = abspath(expanduser(workdir)) if workdir != None else None
  workdir = property( _get_workdir, _set_workdir, 
                      """ Working directory where calculations are performed. 
                      
                          Absolute path is determined when workdir is set.
                      """ )

  def _add_bond(self, args):
    """ Adds/Modifies the bond parameter dictionary.
    
        @param args: - the first element of arg is one endpoint of the bond (say "In"), 
                     - the second element is another endpoint (say "N"),
                     - the last element is a sequence with at most 5 elements
                       defining the bond-stretching parameters (order 2 through 6).
        @type args: "type1", "type2", numpy array
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
    self._add_bond( (A,B, d0, a2, a3, a4, a5, a6) )

  add_bond = add_setter(_add_bond, _add_bond.__doc__)



  def _add_angle(self, args):
    """ Adds/Modifies the angle parameter dictionary.
    
        @param args: - the first element of arg is one endpoint of the angle (say "In"), 
                     - the second element is middle of the angle (say "N"),
                     - the third element is the other endpoint (say "Ga"),
                     - the last element is a sequence with at most 7 elements
                       defining the bond-angle parameters. The first is the
                       gamma parameter, the second the sigma, and the others
                       the bond bending proper.
        @type args: "type1", "type2", numpy array
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
    assert len(params) <= 7 and len(params) > 0, RuntimeError("To few or too many parameters: %s." % (params))
    if name in self.angles: # replaces none value with known value
      for old, new in zip(self.angles[name], params):
        if old == None: old = new
    if str(params[0]).lower()[:3] == "tet": params[0] = -1e0/3e0
    if len(params) < 7: params.extend( 0e0 for u in range(7-len(params)) )

    self.angles[name] = params
    
  def set_angle(self, A, B, gamma = None, sigma = None, a2 = None, a3 = None, a4 = None, a5 = None, a6 = None):
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
    self._add_angle( (A,B, C, gamma, sigma, a2, a3, a4, a5, a6) )

  add_angle = add_setter(_add_angle, _add_angle.__doc__)

  def __repr__(self):
    result  = repr(self.lattice)
    result += repr(self.minimizer)
    result += "# Vff definition.\n"
    result += "vff = %s()\n" % (self.__class__.__name__)
    result += "vff.minimizer = minimizer\n"
    result += "vff.lattice = lattice\n"
    result += "vff.OUTCAR = \"%s\"\n" % (self.OUTCAR)
    result += "vff.ERRCAR = \"%s\"\n" % (self.ERRCAR)
    result += "vff.direction = %s\n" % (repr(self.direction))
    result += "vff.relax = %s\n" % ("True" if self.relax else "False")
    for name, params in self.bonds.items():
      A, B = name.split("-")
      result += "vff.add_bond = \"%s\", \"%s\", (%e" % (A, B, params[0])
      for i, u in enumerate(params[1:]): 
        if sum([abs(j) for j in params[i+1:]]) / float(len(params)-i-1) < 1e-12: break;
        result += ", %e" % (u)
      result += ")\n"
    for name, params in self.angles.items():
      A, B, C = name.split("-")
      result += "vff.add_angle = \"%s\", \"%s\", \"%s\", " % (A, B, C)
      if abs(params[0] + 1e0/3e0) < 1e-12: result += "(\"tet\""
      else: result += "(%e" % (params[0])
      for i, u in enumerate(params[1:]):
        if sum([abs(j) for j in params[i+1:]]) / float(len(params)-i-1) < 1e-12: break;
        result += ", %e" % (u)
      result += ")\n"
    result += "# End of vff definition.\n"
    return result

  def _cout(self, comm):
    """ Creates output name. """
    if self.OUTCAR == None: return "/dev/null"
    return self.OUTCAR if comm.rank == 0 else self.OUTCAR + "." + str(comm.rank)
  def _cerr(self, comm):
    """ Creates error name. """
    if self.ERRCAR == None: return "/dev/null"
    return self.ERRCAR if comm.rank == 0 else self.ERRCAR + "." + str(comm.rank)

  def __call__(self, structure, outdir, comm = None, **kwargs):
    """ Performs calculation """
    import time
    from copy import deepcopy
    from os.path import exists, isdir, abspath, expanduser
    from boost.mpi import world
    from ..opt.changedir import Changedir
    from boost.mpi import world, broadcast

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
    outdir = abspath(expanduser(outdir))

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

    # checks if outdir contains a successful run.
    if broadcast(comm, exists(outdir) if comm.rank == 0 else None, 0):
      extract = Extract(comm = comm, directory = outdir, vff = this)
      if extract.success: return extract # in which case, returns extraction object.
      comm.barrier() # makes sure directory is not created by other proc!
    
    with Changedir(outdir, comm = comm) as current_dir:
      # redirects C/C++/fortran streams
      cout, cerr = this._cout(comm), this._cerr(comm)
      # write stuff to output file
      with open(cout, "w") as file:
        # now for some input variables
        print >> file, "# VFF calculation on ", time.strftime("%m/%d/%y", local_time),\
                       " at ", time.strftime("%I:%M:%S %p", local_time)
        if len(structure.name) != 0: print "# Structure named ", structure.name 
        print >> file, repr(this)
        print >> file, "# Performing VFF calculations. "
        # then calculations
      result, stress = this._run(structure, comm)
      # must close/reopen redirection context, otherwise it seem that they 
      # are closed anyhow on exiting from the Cpp function call. The following context
      with open(cout, "a") as file:
        # finally, the dam results.
        print >> file, "# Result of VFF calculations. "
        print >> file, repr(result)
        print >> file, "stress: (%e, %e, %e),\\\n        (%e, %e, %e),\\\n        (%e, %e, %e)"\
                       % tuple(stress.flat)
        timing = time.time() - timing
        hour = int(float(timing/3600e0))
        minute = int(float((timing - hour*3600)/60e0))
        second = (timing - hour*3600-minute*60)
        print >> file, "# Computed VFF in: %i:%i:%f."  % (hour, minute, second) 


    # checks if result was successful
    extract = Extract(comm = comm, directory = outdir, vff = this)
    assert extract.success, RuntimeError("VFF calculation did not complete in %s.\n" % (outdir))

    return extract


  def _create_functional(self, structure, comm):
    """ Creates the vff functional using cpp extension. """
    from tempfile import NamedTemporaryFile
    from os import remove
    from boost.mpi import broadcast
    from vff import Vff, LayeredVff

    # Creates temporary input file and creates functional
    functional = None
    if comm.rank == 0:
      with NamedTemporaryFile(dir=self.workdir, delete=False) as file: 
        file.write(self._create_input(structure, comm))
      name = broadcast(comm, file.name, root=0)
    else: name = broadcast(comm, root=0) # syncs all procs to make sure we are reading from same file.

    
    functional = Vff(name, comm) if self.direction == None else LayeredVff(name, comm)

    if comm.rank == 0: remove(file.name)

    return functional
     
  def _run(self, structure, comm):
    """ Performs actual calculation. """
    from vff import Vff, LayeredVff
    from ..opt import redirect_all

    # Saves global lattice if set.
    old_lattice = None
    try: old_lattice = structure.lattice
    except RuntimeError: pass
    self.lattice.set_as_crystal_lattice()
    
    cout, cerr = self._cout(comm), self._cerr(comm)
    with open(cerr, "w") as file: pass # file has not yet been opened
    with redirect_all(output=cout, error=cerr, append="True") as oestream:
      print 0
      functional = self._create_functional(structure, comm)
      print 1
      # now performs call
      result, stress = functional(structure, doinit=True)
    
    # unsets lattice.
    if old_lattice != None: old_lattice.set_as_crystal_lattice()
    return result, stress

  def _create_input(self, structure, comm):
    """ Creates a temporary file with input to vff. """

    result = "<?xml version=\"1.0\" standalone=\"no\" ?>\n"\
             "<Job>\n<Functional type=\"vff\""
    if self.direction != None: result += " direction=\"%s %s %s\"" % tuple(self.direction.flat)
    result += ">\n"
    result += "<Minimizer type=\"%s\" tolerance=%e itermax=%i linetolerance=%e\n"\
              "           linestep=%e strategy=\"%s\" verbose=\"%s\" uncertainties=%e\n"\
              "           zeps=%e up=%i gradient=\"%s\"/>\n"\
              % ( self.minimizer.type, self.minimizer.tolerance, self.minimizer.itermax,\
                  self.minimizer.linetolerance, self.minimizer.linestep, self.minimizer.strategy,\
                  "true" if self.minimizer.verbose else "false", self.minimizer.uncertainties,\
                  self.minimizer.uncertainties, self.minimizer.up, \
                  "true" if self.minimizer.use_gradient else "false" )
    types = set(u for site in self.lattice.sites for u in site.type)
    bonds = set(u for name in self.bonds.keys() for u in name.split('-'))
    assert types <= bonds, "Species in structure and vff-input do not match."
    for name, params in self.bonds.items():
      if set(name.split('-')) <= types:
        p = name.split('-')
        p.extend([u for u in params])
        result += "<Bond A=\"%s\" B=\"%s\" d0=%e alpha=%e alpha3=%e alpha4=%e alpha5=%e alpha6=%e />\n"\
                  % tuple(p) 
    for name, params in self.angles.items():
      if set(name.split('-')) <= types:
        p = name.split('-')
        p.extend([u for u in params])
        result += "<Angle A=\"%s\" B=\"%s\" C=\"%s\" gamma=\"%s\" sigma=\"%s\" \n"\
                  "       beta=%e beta3=%e beta4=%e beta5=%e beta6=%e />\n"\
                  % tuple(p) 
    result += "</Functional>\n</Job>"
    return result

