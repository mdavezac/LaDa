""" Valence Force Field functional for zinc-blende materials """

# imports C++ extension
from ..opt.decorators import add_setter
from sys import exit
import sys

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

class Extract(object):
  def __init__(self, directory = "", comm = None, vff = None):
    self.directory = directory
    """ Directory where to check for output. """
    self.comm = comm
    """ Mpi group communicator. """
    self.OUTCAR = vff.OUTCAR
    """ Filename of the OUTCAR file from VASP.
     
        Data will be read from directory/OUTCAR. 
    """
    if vff != None: self.OUTCAR = vff.OUTCAR

    
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
  def success(self):
    """ Checks for VFF success.
        
        At this point, checks for files and 
    """
    from os.path import isdir, exists, join
    import re
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    if not exists(path): return False

    with open(path, "r") as file:
      regex = re.compile(r"""Computed in:""", re.X)
      for line in file:
        if regex.search(line) != None: return True
    return False



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
    self.OUTCAR = "vffout" 
    """ File where to redirect output. """
    self.ERRCAR = "vfferr" 
    """ File where to redirect errors. """
    self.ESCANCAR = "atomic_input"
    """ File where to print atomic inputs for escan. """
    self.SCRIPTCAR = "vff_script.py"
    """ VFF script file for easy reruns. """
    self.workdir = workdir
    """ Working directory. """

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
    if args[0] > args[2]: name = "%s-%s" % (args[2], args[1], args[0])

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

  def __str__(self):
    result  = "%s\n" % (self.lattice)
    result += "# VFF definition:\n"
    result += "vff = Vff()\n"
    result += "vff.lattice = lattice\n"
    result += "vff.OUTCAR = \"%s\"\n" % (self.OUTCAR)
    result += "vff.ERRCAR = \"%s\"\n" % (self.ERRCAR)
    result += "vff.ESCANCAR = \"%s\"\n" % (self.ESCANCAR)
    result += "vff.SCRIPTCAR = \"%s\"\n" % (self.SCRIPTCAR)
    result += "vff.direction = %s\n" % (self.direction)
    result += "vff.relax = %s\n" % ("True" if self.relax else "False")
    for name, params in self.bonds.items():
      A, B = name.split("-")
      result += "vff.add_bond = %s, %s, (%e" % (A, B, params[0])
      for i, u in enumerate(params[1:]): 
        if sum([abs(j) for j in params[i+1:]]) / float(len(params)-i-1) < 1e-12: break;
        result += ", %e" % (u)
      result += ")\n"
    for name, params in self.angles.items():
      A, B, C = name.split("-")
      result += "vff.add_angle = %s, %s, %s, " % (A, B, C)
      if abs(params[0] + 1e0/3e0) < 1e-12: result += "(\"tet\""
      else: result += "(%e" % (params[0])
      for i, u in enumerate(params[1:]):
        if sum([abs(j) for j in params[i+1:]]) / float(len(params)-i-1) < 1e-12: break;
        result += ", %e" % (u)
      result += ")\n"
    result += "%s\n" % (self.minimizer)
    return result


  def __call__(self, structure, outdir, comm = None, **kwargs):
    """ Performs calculation """
    import time
    from copy import deepcopy
    from os.path import exists, isdir, abspath
    from boost.mpi import world
    from ..opt import redirect
    from ..opt.changedir import Changedir

    # bull shit. 
    assert len(self.lattice.sites) == 2, RuntimeError("Lattice is not zinc-blend")
    assert len(self.lattice.sites[0].type) > 0, RuntimeError("No atomic species given in lattice site 0.")
    assert len(self.lattice.sites[1].type) > 0, RuntimeError("No atomic species given in lattice site 1.")
    assert len(self.lattice.sites[1].type) < 3, RuntimeError("Lattice is not binary, pseudo-binary or quaternary.")
    assert len(self.lattice.sites[0].type) < 3, RuntimeError("Lattice is not binary, pseudo-binary or quaternary.")
    assert len(self.lattice.sites[1].type) < 3, RuntimeError("Lattice is not binary, pseudo-binary or quaternary.")
    assert len(self.lattice.sites[0].type) != 1 or  len(self.lattice.sites[1].type) != 2

    timing = time.time() 
    local_time = time.localtime() 

    # gets absolute path.
    outdir = abspath(outdir)

    # make this functor stateless.
    this      = deepcopy(self)
    outdir    = deepcopy(outdir)

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to use for calculations launch. 
    # For added simplicity, minimizer attributes are also directly available.
    for key in kwargs.keys():
      if hasattr(this, key): getattr(this, key).value = kwargs[key]
      elif hasattr(this.minimizer, key): getattr(this.minimizer, key).value = kwargs[key]
      else: raise NameError( "%s attribute unknown of vff or vff.minimizer." % (key) )

    # First checks if directory outdir exists (and is a directory).
    if exists(outdir):
      if not isdir(outdir): raise IOError, "%s exists but is not a directory.\n" % (outdir)
      # checks if it contains a successful run.
      extract = Extract(comm = comm, directory = outdir, vff = this)
      if extract.success: return extract # in which case, returns extraction object.
    
    # Otherwise, performs calculations.
    with Changedir(outdir) as current_dir:
      
      if comm.rank == 0: # a simple script to rerun things.
        with open(this.SCRIPTCAR, "w") as script:
          print >> script, "from boost.mpi import world"
          print >> script, "from lada.vff import Vff"
          print >> script, "from lada.crystal import Structure, Lattice, Atom\n"
          print >> script, this.lattice
          print >> script, "lattice.find_space_group()"
          print >> script, this
          print >> script, structure
          print >> script, "\n# calls vff."
          print >> script, "vff.comm = world"
          print >> script, "vff(structure)"

      # redirects C/C++/fortran streams
      cout = this.OUTCAR if (comm.rank == 0 and len(this.OUTCAR)) else this.OUTCAR + "." + str(comm.rank)
      cerr = this.ERRCAR if (comm.rank == 0 and len(this.ERRCAR)) else this.ERRCAR + "." + str(comm.rank)
      with redirect(cout=cout, cerr=cerr, fout=cout, ferr=cerr) as oestreams:
        # now for some input variables
        print "# VFF calculation on ", time.strftime("%m/%d/%y", local_time),\
              " at ", time.strftime("%I:%M:%S %p", local_time)
        if len(structure.name) != 0: print "# Structure named ", structure.name 
        print this
        print "# Performing calculations. "
        # then calculations
        result, stress = this._run(structure, comm)
        # finally, the dam results.
        print "# Result of vff calculations. "
        print result
        print "stress: (%e, %e, %e),\\\n        (%e, %e, %e),\\\n        (%e, %e, %e)"\
              % tuple(stress.flat)
        timing = time.time() - timing
        hour = timing/3600
        minute = (timing - hour*3600)/60
        second = (timing - hour*3600-minute*60)
        print "# Computed in: %i:%i:%f."  % (int(hour), int(minute), second) 

    # checks if result was successful
    extract = Extract(comm = comm, directory = outdir, vff = this)
    assert extract.success, RuntimeError("VFF calculation did not complete in %s.\n" % (outdir))

    return extract


  def _run(self, structure, comm):
    """ Performs actual calculation. """
    from vff import Vff, LayeredVff
    from tempfile import NamedTemporaryFile
    from os import remove

    # Saves global lattice if set.
    old_lattice = None
    try: old_lattice = structure.lattice
    except RuntimeError: pass
    self.lattice.set_as_crystal_lattice()
    
    # Creates temporary input file and creates functional
    functional = None
    if comm.rank == 0:
      with NamedTemporaryFile(dir=self.workdir, delete=False) as file: 
        file.write(self._create_input(structure, comm))
    comm.barrier() # syncs all procs to make sure we are reading from same file.
    functional = Vff(file.name, comm) if self.direction == None else LayeredVff(file.name, comm)
    if comm.rank == 0: remove(file.name)


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
    types = set(u.type for u in structure.atoms)
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

