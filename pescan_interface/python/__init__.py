""" Interface module for pescan. """
from ..opt import __load_pescan_in_global_namespace__
if __load_pescan_in_global_namespace__:
  from DLFCN import RTLD_NOW as _RTLD_NOW, RTLD_GLOBAL as _RTLD_GLOBAL
  from sys import getdlopenflags as _getdlopenflags, setdlopenflags as _setdlopenflags
  flags = _getdlopenflags()
  _setdlopenflags(_RTLD_NOW|_RTLD_GLOBAL)
  import _escan
  _setdlopenflags(flags)
else: import _escan
from bandstructure import band_structure

def call_escan(comm, atom="atom_input", pot="pot_input", escan="escan_input"):
  """ Calls escan functional in current directory.

      Before calling the functional, the files are propagated such that each
      proc can read its own. What is more, the inputs are modified such that
      L{pot} does refer to L{atom}, and L{escan} to both atomic config and
      potential files.
      @param comm: Processes on which to execute.
      @type comm: boost.mpi.Communicator or mpi4py
      @param atom: file with atomic input (from vff). 
      @param pot: file with input to potential generation.
      @param escan: file with input to escan itself.
  """
  from os import remove
  from boost.mpi import world, broadcast
  from _escan import _call_escan, _call_genpot
  
  # escan is free for all, with every proc reading at the same time.
  # hence it must be from different files.
  atominput = "%s.%i" % (atom, world.rank)
  potinput = "%s.%i" % ("pot.input", world.rank)
  potoutput = "%s.%i" % ("pot.output", world.rank)
  escaninput = "%s.%i" % ("escan_input", world.rank)

  # propagates atomic input first
  buffer = None
  if comm.rank == 0:
    with open(atom, "r") as file: buffer = file.read()
  buffer = broadcast(comm, buffer, 0)
  with open(atominput, "w") as file: file.write(buffer)

  # propagates pot input second
  buffer = None
  if comm.rank == 0:
    with open(pot, "r") as file: buffer = file.readlines()
    buffer[0] = buffer[0].replace( buffer[0].split()[0], atom )
  buffer = broadcast(comm, buffer, 0)
  with open(potinput, "w") as file: file.writelines(buffer)

  # propagates pot input second
  buffer = None
  if comm.rank == 0:
    with open(escan, "r") as file: buffer = file.readlines()
    buffer[0] = buffer[0].replace( buffer[0].split()[1], potoutput )
    buffer[12] = buffer[12].replace( buffer[12].split()[1], atominput )
  buffer = broadcast(comm, buffer, 0)
  with open(escaninput, "w") as file: file.writelines(buffer)

  #  calls escan at last. 
  _call_genpot(comm)
  _call_escan(comm)

  remove(atominput)
  remove(potinput)
  remove(escaninput)

folded_spectrum = 0
""" Folded spectrum method. """
all_electron = 1
""" All electron method. """
local_hamiltonian = 0
""" Local hamiltonian. """
nonlocal_hamiltonian = 1
""" Non-local hamiltonian. """
spin_orbit_hamiltonian = 2
""" Spin-orbit hamiltonian. """

def AtomicPotentials(object):
  """ Holds parameters to atomic potentials. """
  def __init__( self, filepath, specie, nonlocal=None, s=None, p=None,\
                d=None, pnl=None, dnl=None, Z = None ):
    import physics
    from os.path import abspath

    self.filepath = abspath(filepath)
    self.specie = specie
    self.nonlocal = None if nonlocal == None else abspath(nonlocal)
    self.s =  s
    """ s parameter """
    self.p =  p
    """ p parameter """
    self.d =  d
    """ d parameter """
    self.pnl = pnl
    """ p non-local parameter """
    self.dnl = dnl
    """ d non-local parameter """
    self.Z = physics.Z(self.specie) if Z == None else Z
    """ Atomic number of the specie. 

        If not given when initialized, will use atomic specie database.
    """

class Pescan(object):
  """ Performs PESCAN calculations, from structure relaxation to wavefunctions. """

  def __init__(self, workdir):
    """ Initializes a PESCAN functional. """
    from numpy import zeros
    from ..vff import Vff

    super(Pescan, self).__init__()
    self.vff = Vff() 
    """ The L{Vff} functional with which to relax a structure. """
    self.OUTCAR = "escan_out" 
    """ Escan output file. """
    self.ERRCAR = "escan_err"
    """ Escan error file. """
    self.WAVECAR = "wavefunctions"
    """ Wavefunction file (in g-space). """
    self.eref = None
    """ Reference energy for folded spectrum method.
    
        Set to None for all electron diagonalization.
    """
    self.cutoff = 8.2
    """ Cutoff energy for plane-wave expansion. """
    self.smooth = 1
    """ Smooth potential scaling. """
    self.kinetic_scaling = 1
    """ Smooth kinetic energy scaling. """
    self.nbstates = 1
    """ Number of states to compute. """
    self.itermax = 20
    """ Maximum number of iterations when minimzing wavefuntions. """
    self.nlines = 50
    """ Conjugate gradient line minimization parameter. """
    self.tolerance = 1e-12
    """ Requested minimization accuracy. """
    self.input_wavefunctions = []
    """ Indices of wavefunctions to read in. """
    self.INWAVECAR = "input_escan_wavefunctions"
    """ Filename of input wavefunctions. """
    self.kpoint = zeros((3,1), dtype="float64")
    """ k-point at which to perform calculations. """
    self.potential = spin_orbit_hamiltonian
    """ Type of hamiltonian to use. """
    self.rspace_cutoff = 5
    """ real-space projector cutoff. """
    self.atomic_potentials = None
    """ Parameters to atomic potentials. """
    self.fft_mesh = (18, 18, 18)
    """ Fourrier Transform mesh. """
    self.dnc_mesh = None
    """ Divide and conquer mesh.
     
        Set to None if no divide and conquer is required. 
    """
    self.overlap_mesh = None
    """ Overlap of divide and conquer mesh. 
     
        Set to None if no divide and conquer is required. 
    """
    self.do_genpot = True
    """ If true, the potential is generated.
          
        Otherwise, the file L{_POTCAR} from the output directory is used. 
    """
    self.do_relax = True
    """ If true, the structure is relaxed using vff.
          
        Otherwise, the file L{_POSCAR} from the output directory is used. 
    """

    self._ATOMPINPUT = "atomic_input"
    """ Private reference to the atomic input file. """
    self._workdir = workdir
    """ Private reference to the working directory. """
    self._POTCAR = "potential"
    """ Private reference to the potential file generated by genpot. """
    self._WAVECAR = "gspace_wavefunctions"
    """ Private reference to the wavefunctions file generated by escan. """
    self._MASKCAR = "maskr"
    """ Private reference to the maskr projector file. """
    self._INCAR = "escan.input"
    """ Private reference to the maskr projector file. """


  def _get_workdir(self): return self._workdir
  def _set_workdir(self, workdir):
    from os.path import abspath
    self._workdir = abspath(workdir) if workdir != None else None
  workdir = property( _get_workdir, _set_workdir, 
                      """ Working directory where calculations are performed. 
                      
                          Absolute path is determined when workdir is set.
                      """ )

  def _get_maskr(self): return self._maskr
  def _set_maskr(self, workdir): 
    from os.path import abspath
    self._maskr = abspath(workdir) if workdir != None else None
  maskr = property( _get_maskr, _set_maskr, 
                      """ maskr file for projectors.
                      
                          Absolute path is determined when maskr is set.
                      """
                    )

  def _get_atompots(self): return self._atompots
  def _set_atompots(self): 
    from os.path import abspath
    self._maskr = abspath(workdir) if workdir != None else None
  maskr = property( _get_atompots, _set_atompots, 
                    """ Files to the atomic potentials.
                    
                        Absolute path is determined when the atomic potentials are set.
                    """
                   )

  def _add_potential(self, args):
    """ Adds an atomic potential to escan. """
    assert len(args) > 2, RuntimeError("Atomic  potentials need at least two parameters.")
    assert len(args) < 9, RuntimeError("Too many parameters when setting atomic potentials.")
    result = AtomicPotentials(*args)
    self.atomic_potentials( AtomicPotentials(*args) )
  add_potential = add_setter( _add_potential,\
                              """ Adds atomic potential to escan.
                                  
                                  - first argument is the path to the atomic potential. 
                                    The absolute path is deduced when set.
                                  - second argument is the symbol of the atomic specie.
                                  - third argument is the path to the non-local potential file. 
                                    If None, then no non-local argument is added. Defaults to None.
                                  - fourth trough eigths arguments are the
                                    non-local potential parameters s, p, d,
                                    pnl, dnl. Defaults to None (eg 0).
                                  - last argument is atomic number. Defaults to L{physics.Z}(specie).
                              """ )

  def __str__(self):
    result  = str(self.vff)
    result += "# Pescan definition."
    result += "escan.eref                  = %s\n" % ( "None" if self.eref == None else str(self.eref) )
    result += "escan.cutoff                = %f\n" % (self.cutoff)
    result += "escan.smooth                = %f\n" % (self.smooth)
    result += "escan.kinetic_scaling       = %f\n" % (self.kinetic_scaling)
    result += "escan.nbstates              = %i\n" % (self.nbstates)
    result += "escan.itermax               = %i\n" % (self.itermax)
    result += "escan.nlines                = %i\n" % (self.nlines)
    result += "escan.tolerance             = %e\n" % (self.tolerance)
    result += "escan.rspace_cutoff         = %f\n" % (self.rspace_cutoff)
    result += "escan.fft_mesh              = %i, %i, %i\n" % self.fft_mesh
    result += "escan.do_genpot             = %s" % ("True" if self.do_genpot else "False")
    result += "escan.do_relax              = %s" % ("True" if self.do_relax else "False")
    result += "escan.input_wavefunctions = " % (rep(self.input_wavefunctions))
    result += "escan.INWAVECAR = " % (self.INWAVECAR)
    result += "escan.kpoint                = " % (rep(self.kpoint))
    result += "escan.dnc_mesh              = " % (rep(self.dnc_mesh))
    result += "escan.overlap_mesh          = " % (rep(self.overlap_mesh))

    result += "escan.potential             = " % (self.potential)
    result += "escan.atomic_potentials     = " % (self.atomic_potentials)


    self.OUTCAR = "escan_out" 
    self.ERRCAR = "escan_err"
    self.WAVECAR = "wavefunctions"

    self._ATOMPINPUT = "atomic_input"
    """ Private reference to the atomic input file. """
    self._workdir = workdir
    """ Private reference to the working directory. """
    self._POTCAR = "potential"
    """ Private reference to the potential file generated by genpot. """
    self._maskr = "maskr"
    """ Private reference to the maskr projector file. """
    self._GENCAR = "pot.input"


  def __call__(self, structure, outdir = None, workdir = None,  comm = None, **kwargs):
    """ Performs calculation """
    import time
    from copy import deepcopy
    from os.path import exists, isdir, abspath
    from boost.mpi import world
    from ..opt import redirect_all, redirect
    from ..opt.changedir import Changedir

    timing = time.time() 
    local_time = time.localtime() 

    # make this functor stateless.
    this      = deepcopy(self)
    outdir    = deepcopy(outdir)

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to use for calculations launch. 
    for key in kwargs.keys():
      if hasattr(this, key): getattr(this, key).value = kwargs[key]
      else: raise NameError( "%s attribute unknown of vff or vff.minimizer." % (key) )

    # First checks if directory outdir exists (and is a directory).
    if exists(outdir):
      if not isdir(outdir): raise IOError, "%s exists but is not a directory.\n" % (outdir)
      # checks if it contains a successful run.
      extract = Extract(comm = comm, directory = outdir, vff = this)
      if extract.success: return extract # in which case, returns extraction object.

  def _run_genpot(self, comm, types):
    """ Runs genpot only """
    from ._escan import _call_genpot
    from shutil import copyfile
    from os.path import basename

    assert self.atomic_potentials != None, RuntimeError("Atomic potentials are not set.")
    # Creates temporary input file and creates functional
    with Changedir(self.workdir) as workdir:
      dnc_mesh = self.dnc_mesh if self.dnc_mesh != None else self.mesh
      overlap_mesh = self.overlap_mesh if self.overlap_mesh != None else (0,0,0)
      with open(self._GENCAR + "." + str(comm.rank), "w") as file:
        file.write( "%s\n%i %i %i\n%i %i %i\n%i %i %i"\
                    % (self._POSCAR, *self.mesh, *dnc_mesh, *overlap_mesh) )
        for u in [v in self.atomic_potentials if v.specie in types]:
          # adds to list of potentials
          filepath = basename(u.filepath)
          file.write(filepath + "\n") 
          # copy potential files as well.
          copyfile(u.filepath, filepath)

      comm.barrier() # syncs all procs to make sure we are reading from same file.
      _call_genpot(comm)

  def _run_escan(self, comm, types, scale):
    """ Runs escan only """
    from shutil import copyfile
    from os.path import basename
    from ._escan import _call_pescan

    assert self.atomic_potentials != None, RuntimeError("Atomic potentials are not set.")
    # Creates temporary input file and creates functional
    with Changedir(self.workdir) as workdir:
      with open(self._INCAR + "." + str(comm.rank), "w") as file:
        print >> file, "1 %s.%i" % (self._POTCAR, comm.rank) )
        print >> file, "2 %si" % (self._WAVECAR) )
        print >> file, "3 %i # %s" % ((1, "folded spectrum") if self.eref != None else (2, "all electron"))
        print >> file, "4 %f %f %f %f # Eref, cutoff, smooth, kinetic scaling"\
                       % ( self.eref if self.eref != None else 0,\
                           self.cutoff, self.smooth, self.kinetic )
        print >> file, "5 %i # number of states" % (self.nbstates)
        print >> file, "6 %i %i %e # niter, nline, tolerance" % (self.niter, self.nlines, self.tolerance)

        nowfns = self.input_wavefunctions == None
        if not nowfns: nowfns = len(self.input_wavefunctions) == 0
        if not nowfns: print >> file "7 0 # no input wfns\n8 0 # wfns indices"
        else:
          print >> file "7 %i\n8 %i" % (len(self.input_wavefunctions), self.input_wavefunctions[0])
          for u in self.input_wavefunctions[1:]:
            print >> file, str(u),
          print >> file
        print >> file, "9 %s # input wavefunction filename" % (self.INWAVECAR)

        print >> file, "10 0 1 1 1 whatever"

        if np.linalg.norm(self.kpoint) < 1e-12: print >> file, "11 0 0 0 0 0"
        else: print >> file, "11 1 %i %i %i" % tuple((self.kpoint * scale / a0("A")).flat)
        
        if   self.potential == local_hamiltonian: print >> file, "12 1 # local hamiltonian" 
        elif self.potential == nonlocal_hamiltonian: print >> file, "12 1 # non-local hamiltonian" 
        elif self.potential == spin_orbit_hamiltonian: print >> file, "12 1 # spin orbit hamiltonian" 
        else: raise RuntimeError("Unknown potential requested.")
        
        print >> file, "13 ", self._POSCAR  
        print >> file, "14 ", self.rspace_cutoff, "# real-space cutoff" 

        soatoms = [u for u in self.atomic_potentials if u in types]\
                  if self.potential == spin_orbit_hamiltonian \
                  else []
        print >> file, "15 ", len(soatoms), "# Number of spin-orbit potentials"
        for i, u in enumerate(soatoms):
          filepath = basename(u.nonlocal)
          print >> file, i + 15, filepath, u.Z, u.s, u.p. u.d, u.pnl, u.dnl
          copyfile(u.nonlocal, filepath)

      comm.barrier() # syncs all procs to make sure we are reading from same file.
      copyfile(self._MASKCAR, basename(self._MASCAR))
      _call_pescan(comm)
