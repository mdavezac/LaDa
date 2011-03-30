""" ESCAN functional wrapper. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Functional', 'folded_spectrum', 'all_electron']
from .. import lada_with_mpi
from ..opt.decorators import add_setter
from ._extract import Extract

folded_spectrum = 0
""" Folded spectrum method. """
all_electron = 1
""" All electron method. """


class Escan(object):
  """ Performs ESCAN calculations, from structure relaxation to wavefunctions. """

  Extract = staticmethod(Extract)
  """ Class for output extraction. """

  def __init__(self, inplace=True, workdir=None):
    """ Initializes a ESCAN functional. """
    from numpy import zeros
    from ..vff import Vff
    from ..opt import RelativeDirectory
    from ._potential import soH

    object.__init__(self)
    self.inplace = inplace
    """ If True calculations are performed in the output directory. """
    # checks inplace vs workdir
    if self.inplace: 
      assert workdir == None, ValueError("Cannot use both workdir and inplace attributes.")

    self.vff = Vff() 
    """ The `lada.vff.Vff` functional with which to relax a structure. """
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
    """ k-point at which to perform calculations.
    
        By default, k-points are given in cartesian coordinates, and units of
        2pi/structure.scale. However, since relaxation will deform the BZ, the
        kpoint given to escan is the one deformed to the relaxed BZ. E.g, the
        meaning of [1,0,0] stays [1,0,0], despite possible relaxation. To turn
        this behavior off, use `do_relax_kpoint`.
    """
    self.potential = soH
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
    self.vffrun = None
    """ If None, the structure is relaxed using vff.
          
        Otherwise, it should be an extraction object returned from a previous
        run where vff was computed.
    """
    self.genpotrun = None
    """ If None, the potential is generated.
          
        Otherwise, it should be the extraction object returned by a previous
        run which computed the potential.
    """
    self.do_escan = True
    """ If true, calculations are performed. """

    self._POSCAR = "atomic_input"
    """ Private reference to the atomic input file. """
    self._workdir = RelativeDirectory(path=workdir)
    """ Working directory where calculations are performed. """
    self._POTCAR = "pot.output"
    """ Private reference to the potential file generated by genpot. """
    self._maskr = RelativeDirectory("maskr")
    """ Reference to the maskr projector file. """
    self._INCAR = "escan_input"
    """ Private reference to the escan input file. """
    self._GENCAR = "pot.input"
    """ Private reference to the genpot input file. """
    self._FUNCCAR = "ESCANCAR"
    """ Private reference to the functional pickle. """
    self.do_relax_kpoint = False
    """ Whether to deform kpoints from input cell to relaxed cell.

        Default is True. Relaxed cell is taken from `_POSCAR`
    """
    self.print_from_all = False
    """ If True, each node will print. """
    self.symlink = True
    """ Prefer symlink to actual copy where possible. """

  @property
  def relax(self):
    """ Whether or not to relax the structure with vff. """
    return self.vff.relax
  @relax.setter
  def relax(self, value): self.vff.relax = value

  @property
  def is_krammer(self):
    """ True if wavefunction is a spinor. """
    from numpy.linalg import norm
    from . import soH
    return True if (norm(self.kpoint) < 1e-12 and self.potential == soH) else False

  @property
  def _dont_deform_kpoint(self):
    from warnings import warn
    warn( DeprecationWarning('_dont_deform_kpoint is deprecated in favor of do_relax_kpoint.' ),
          stacklevel=2 )
    return not self.do_relax_kpoint
  @_dont_deform_kpoint.setter
  def _dont_deform_kpoint(self, value):
    from warnings import warn
    warn( DeprecationWarning('_dont_deform_kpoint is deprecated in favor of do_relax_kpoint.' ),
          stacklevel=2 )
    self.do_relax_kpoint = value

  @property
  def maskr(self): 
    """ Reference to the maskr projector file. """
    return self._maskr.path
  @maskr.setter
  def maskr(self, value): self._maskr.path = value

  @property
  def workdir(self): 
    """ Directory where calculations are performed. 
    
        This variable is used only if inplace is False. By default,
        calculations are inplace.
    """
    return self._workdir.path
  @maskr.setter
  def workdir(self, value): self._workdir.path = value
    
  @property
  def lattice(self):
    """ Lattice to use with escan and vff functionals. """
    return self.vff.lattice
  @lattice.setter
  def lattice(self, value): self.vff.lattice = value

  @add_setter
  def add_potential(self, args):
    """ Adds atomic potential to escan.
        
        This property can only be set and never gotten. It serves to add
        potentials to the escan functional.

        >>> escan.add_potential = *args

        Where *args* is a tuple of variable length containing the following:

        - first argument is the path to the atomic potential. 
          The absolute path is deduced when set.
        - second argument is the path to the non-local potential file. 
          If None, then no non-local argument is added. Defaults to None.
        - third trough seventh arguments are the
          non-local potential parameters s, p, d,
          pnl, dnl. Defaults to None (eg 0).
    """ 
    from ._potential import AtomicPotential
    assert len(args) > 2, RuntimeError("Atomic  potentials need at least two parameters.")
    assert len(args) < 9, RuntimeError("Too many parameters when setting atomic potentials.")
    if self.atomic_potentials == None: self.atomic_potentials = []
    self.atomic_potentials.append( AtomicPotential(*args) )

  def __repr__(self):
    from os.path import relpath
    from ._potential import localH, nonlocalH, soH
    result  = str(self.vff).replace("functional", "vff_functional")
    result += "# Escan definition.\n"
    result += "functional = %s()\n" % (self.__class__.__name__)

    # create format string for public data members.
    max_length, string, _string, values = len('add_potential'), '', '', {}
    for key, value in self.__dict__.items():
      if key[0] == '_': continue
      if key == 'potential': continue
      if key == 'workdir': continue
      if key == 'atomic_potentials': continue
      if key == 'vff': continue
      if key == 'lattice': continue
      try: r = repr(value).rstrip().lstrip()
      except: continue
      if r[0] == '<' or r[-1] == '>': continue
      max_length = max(max_length, len('{0}'.format(key)))
      string += 'functional.{{{0}: <{{_mxlgth_repr_}}}} = {1}\n'.format(key, r)
      values[key] = key
    # create format string for private data members.
    for key, value in self.__dict__.items():
      if key[0] != '_': continue
      if key == '_workdir': continue
      if key == '_maskr': continue
      try: r = repr(value).rstrip().lstrip()
      except: continue
      if r[0] == '<' or r[-1] == '>': continue
      max_length = max(max_length, len('{0}'.format(key)))
      _string += 'functional.{{{0}: <{{_mxlgth_repr_}}}} = {1}\n'.format(key, r)
      values[key] = key

    result += "functional.{0: <{1}} = vff_functional\n".format('vff', max_length)
    values['_mxlgth_repr_'] = max_length
    if self.potential == localH:
      result += "functional.{0: <{1}} = localH\n".format('potential', max_length)
    elif self.potential == nonlocalH:
      result += "functional.{0: <{1}} = nonlocalH\n".format('potential', max_length)
    elif self.potential == soH:
      result += "functional.{0: <{1}} = soH\n".format('potential', max_length)
    else: raise RuntimeError("unknown hamiltonnian {0}.".format(soH))
    for pot in self.atomic_potentials:
      result += "functional.{0: <{1}} = {2}\n".format('add_potential', max_length, repr(pot))
    if self.inplace == False: 
      result += "functional.{0: <{1}} = {2}\n"\
                .format('workdir', max_length, repr(self._workdir.unexpanded))
    result += string.format(**values)
    result += 'functional.{0: <{1}} = {2}\n'.format('maskr', max_length, self._maskr.repr())
    result += _string.format(**values)
    result += "# End of escan definition."

    module = self.__class__.__module__ 
    classname = self.__class__.__name__ 
    header = "from numpy import array\n"\
             "from lada.escan import {0}, soH, localH, nonlocalH\n".format(classname)
    return header + result

  def __call__(self, structure, outdir = None, comm = None, overwrite=False, \
               norun=False, workdir=None, do_vff=True, do_genpot=True, **kwargs):
    """ Performs calculation """
    from copy import deepcopy
    from os import getcwd
    from os.path import exists, isdir, abspath, basename, join, expanduser
    from shutil import rmtree
    from cPickle import dump
    from ..opt import copyfile
    from ..opt.changedir import Changedir
    from ..opt.tempdir import Tempdir
    from ..mpi import Communicator

    comm = Communicator(comm, with_world=True)

    if outdir == None: outdir = getcwd()

    # make this functor stateless.
    this      = deepcopy(self)
    outdir    = abspath(expanduser(outdir))

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to use for calculations launch. 
    # If an attribute cannot be found to exist in escan, then vff attributes
    # are checked, and lastly vff.minimizer attributes.
    for key in kwargs.keys():
      if hasattr(this, key): setattr(this, key, kwargs[key])
      elif hasattr(this.vff, key): setattr(this.vff, key, kwargs[key])
      elif hasattr(this.vff.minimizer, key): setattr(this.vff.minimizer, key, kwargs[key])
      else: raise NameError( "%s attribute unknown of escan." % (key) )

    # checks if outdir contains a successful run.
    does_exist, overwrite = comm.broadcast((exists(outdir) if comm.is_root else None, overwrite))
    if does_exist and not overwrite:
      extract = Escan.Extract(directory=outdir, escan=this, comm=comm)
      if extract.success: return extract # in which case, returns extraction object.
    comm.barrier() # makes sure directory is not created by other proc!

    # changes to temporary working directory
    if self.inplace: context = Changedir(outdir, comm=comm) 
    else:            context = Tempdir(workdir = this.workdir, comm=comm)
    with context as this._tempdir: 

      # Saves FUNCCAR.
      if comm.is_root:
        path = join(abspath(this._tempdir), this._FUNCCAR)
        with open(path, "w") as file: dump(this, file)
  
      # performs calculation.
      this._run(structure, outdir, comm, overwrite, norun, do_vff, do_genpot)
  
      # copies output files.
      if not self.inplace:
        with Changedir(outdir, comm = comm) as cwd:
          for file in  [ this._POSCAR, 
                         this._POTCAR, 
                         this.FUNCCAR, 
                         this._cout(comm), 
                         this._cerr(comm), 
                         this.vff._cout(comm),
                         this.vff._cerr(comm),
                         this.WAVECAR if comm.rank == 0  else None ]:
            copyfile(file, this._tempdir, 'same exists null', None, aslink=True)
  
    return self.Extract(comm = comm, directory = outdir, escan = this)

  def _cout(self, comm):
    """ Creates output name. """
    if self.OUTCAR == None: return "/dev/null"
    if comm.is_root: return self.OUTCAR
    return self.OUTCAR + "." + str(comm.rank) if self.print_from_all else "/dev/null"

  def _cerr(self, comm):
    """ Creates error name. """
    if self.ERRCAR == None: return "/dev/null"
    if comm.is_root: return self.ERRCAR
    return self.ERRCAR + "." + str(comm.rank) if self.print_from_all else "/dev/null"


  def _run(self, structure, outdir, comm, overwrite, norun, do_vff, do_genpot):
    """ Performs escan calculation. """
    import time
    from os.path import join
    from ..opt.changedir import Changedir

    if self.genpotrun != None and self.vffrun != None and self.do_escan == False:
      print "Nothing to do? no relaxation, no genpot, no escan?" 
      return None
    timing = time.time() 
    local_time = time.localtime() 

    # prints some output first
    cout, cerr = self._cout(comm), self._cerr(comm)
    with Changedir(self._tempdir, comm=comm) as cwd:
      with open(cout, "w") as file: 
        print >>file, "# Escan calculation on ", time.strftime("%m/%d/%y", local_time),\
                      " at ", time.strftime("%I:%M:%S %p", local_time)
        if comm.is_mpi:
          from ..mpi import world
          file.write("# Computing with {0} processors of {1}.\n".format(comm.size, world.size))
        if len(structure.name) != 0: file.write("# Structure named {0}.".format(structure.name))
        # changes directory to get relative paths.
        with Changedir(outdir, comm=comm) as outdir_wd:
          print >>file, repr(self)
        print >>file, "# Performing calculations. "
      
      # makes calls to run
      if do_vff: self._run_vff(structure, outdir, comm, cout, overwrite, norun)
      local_comm = self._local_comm(comm)
      if local_comm != None:
        if do_genpot: self._run_genpot(local_comm, outdir, norun)
        if self.do_escan: self._run_escan(local_comm, structure, norun)

      # don't print timeing if not running.
      if norun == True: return


      with open(cout, "a") as file: 
        timing = time.time() - timing
        hour = int(float(timing/3600e0))
        minute = int(float((timing - hour*3600)/60e0))
        second = (timing - hour*3600-minute*60)
        file.write("# Computed ESCAN in: %i:%i:%f.\n"  % (hour, minute, second))
      
      if self.do_escan: 
        extract = Extract(comm=comm, directory = outdir, escan = self)
        assert extract.success, RuntimeError("Escan calculations did not complete.")

  def _local_comm(self, comm):
    """ Communicator over which calculations are done. """
    if not comm.is_mpi: return comm
    fftsize = self.fft_mesh[0] * self.fft_mesh[1] * self.fft_mesh[2]
    for m in range(comm.size, 0, -1):
      if fftsize % m == 0: break
    norun = comm.rank >= m
    result = comm.split(0 if norun else 1)
    return None if norun else result

  def _run_vff(self, structure, outdir, comm, cout, overwrite, norun):
    """ Gets atomic input ready, with or without relaxation. """
    from os.path import join, samefile, exists
    from ..vff import Extract as ExtractVff
    from ..opt import copyfile

    if comm.is_root and self.vffrun != None:
      vffrun = self.vffrun.solo()
      POSCAR = join(vffrun.directory, vffrun.functional._POSCAR)
      rstr = vffrun.structure
      if exists(POSCAR): copyfile(POSCAR, self._POSCAR, 'same', None, self.symlink)
      else: vffrun.write_escan_input(self._POSCAR, rstr)
      VFFCOUT = vffrun.functional.vff._cout(comm)
      VFFCOUT = join(vffrun.directory, VFFCOUT)
      copyfile(VFFCOUT, self.vff._cout(comm), 'same exists null', None, self.symlink)

    if self.vffrun != None or norun == True: return
    
    out = self.vff(structure, outdir=outdir, comm=comm, overwrite=overwrite)
    assert out.success, RuntimeError("VFF relaxation did not succeed.")
    if comm.is_root: out.solo().write_escan_input(self._POSCAR, out.solo().structure)

    # copies vff output to stdout. This way, only one outcar.
    if comm.is_root and out.OUTCAR != self.OUTCAR:
      s = out.solo()
      with open(join(s.directory, s.OUTCAR)) as file_in: 
        with open(cout, "aw") as file_out: 
          for line in file_in:
            if line.find("# VFF calculation on ") != -1: print >>file_out, line[:-1]
            if line == "# Performing VFF calculations. ": break
          print >>file_out, line[:-1]
          for line in file_in:
            if line.find("# Computed VFF in:") != -1: break
            print >>file_out, line[:-1]
          print >>file_out, line[:-1]


  def _run_genpot(self, comm, outdir, norun):
    """ Runs genpot only """
    from os.path import basename, exists, join, samefile
    from ..opt import redirect, copyfile

    # using genpot from previous run
    if comm.is_root and self.genpotrun != None:
      genpotrun = self.genpotrun.solo()
      POTCAR = join(genpotrun.directory, genpotrun.functional._POTCAR)
      potcar = self._POTCAR
      copyfile(POTCAR, potcar, 'same exists', None, self.symlink)
      copyfile(self.maskr, nothrow='same', comm=None, symlink=self.symlink)
      for pot in self.atomic_potentials:
        copyfile(pot.nonlocal, nothrow='none same', comm=None, symlink=self.symlink)
    if self.genpotrun != None: return

    assert self.atomic_potentials != None, RuntimeError("Atomic potentials are not set.")
    # Creates temporary input file and creates functional
    dnc_mesh = self.dnc_mesh if self.dnc_mesh != None else self.fft_mesh
    overlap_mesh = self.overlap_mesh if self.overlap_mesh != None else (0,0,0)
    if comm.is_root: 
      with open(self._GENCAR, "w") as file:
        file.write( "%s\n%i %i %i\n%i %i %i\n%i %i %i\n%f\n%i\n"\
                    % ( self._POSCAR, self.fft_mesh[0], self.fft_mesh[1], self.fft_mesh[2], \
                        dnc_mesh[0], dnc_mesh[1], dnc_mesh[2],\
                        overlap_mesh[0], overlap_mesh[1], overlap_mesh[2], self.cutoff,\
                        len(self.atomic_potentials) ))
        for pot in self.atomic_potentials:
          # adds to list of potentials
          file.write(basename(pot.filepath) + "\n") 
    for pot in self.atomic_potentials:
      # copy potential files as well.
      copyfile(pot.filepath, nothrow='same', comm=comm, symlink=self.symlink)
      copyfile(pot.nonlocal, nothrow='same None', comm=comm, symlink=self.symlink)

    copyfile(self.maskr, nothrow='same', comm=comm, symlink=self.symlink)

    if norun == False:
      with redirect(fout=self._cout(comm), ferr=self._cerr(comm), append=True) as oestreams: 
        assert comm.real, RuntimeError('Cannot run escan without mpi.')
        from ._escan import _call_genpot
        _call_genpot(comm)


  def _write_incar(self, comm, structure, norun=False):
    """ Writes escan input to file. """
    from os.path import basename
    from numpy.linalg import norm
    from quantities import eV
    from ..physics import Ry
    from ._potential import soH, nonlocalH, localH
    assert self.atomic_potentials != None, RuntimeError("Atomic potentials are not set.")
    # Creates temporary input file and creates functional
    kpoint = (0,0,0,0,0) if norm(self.kpoint) < 1e-12\
             else self._get_kpoint(structure, comm, norun)
    if comm.is_root:
      with open(self._INCAR, "w") as file:
        file.write('1 {0}\n'.format(self._POTCAR))
        file.write('2 {0.WAVECAR}\n3 {1}\n'.format(self, 1 if self.eref != None else 2) )
        eref = self.eref if self.eref != None else 0
        if hasattr(eref, "rescale"): eref = float(eref.rescale(eV))
        cutoff = self.cutoff
        if hasattr(cutoff, "rescale"): cutoff = float(cutoff.rescale(Ry))
        file.write('4 {0} {1} {2.smooth} {2.kinetic_scaling}\n'.format(eref, cutoff, self))
        nbstates = self.nbstates
        if self.potential != soH or norm(self.kpoint) < 1e-6: nbstates = max(1, self.nbstates/2)
        assert nbstates > 0, ValueError("Cannot have less than 1 state ({0}).".format(nbstates))
        file.write( '5 {0}\n6 {1.itermax} {1.nlines} {1.tolerance}\n'.format(nbstates, self))
        nowfns = self.input_wavefunctions == None
        if not nowfns: nowfns = len(self.input_wavefunctions) == 0
        if nowfns: file.write('7 0\n8 0\n')
        else:
          file.write( '7 {0}\n8 {1}'\
                      .format(len(self.input_wavefunctions), self.input_wavefunctions[0]) )
          for u in self.input_wavefunctions[1:]: file.write(' {0}'.format(u))
          file.write('\n')
        file.write('9 {0.INWAVECAR}\n10 0 1 1 1 0\n11 {1[0]} {1[1]} {1[2]} {1[3]} {1[4]}\n'\
                   .format(self, kpoint))
        
        if   self.potential == localH:    file.write("12 1 # local hamiltonian\n")
        elif self.potential == nonlocalH: file.write("12 2 # non-local hamiltonian\n")
        elif self.potential == soH:       file.write("12 3 # spin orbit hamiltonian\n")
        else: raise RuntimeError("Unknown potential requested.")
        
        file.write('13 {0}\n'.format(self._POSCAR))
        file.write('14 {0.rspace_cutoff}\n'.format(self))
  
        if self.potential != soH: file.write('15 0\n')
        else:
          file.write('15 {0}\n'.format(len(self.atomic_potentials)))
          for i, pot in enumerate(self.atomic_potentials):
            filepath = basename(pot.nonlocal)
            file.write( '{0} {1} {2} {3} {4} {5} {6} {7}\n'\
                        .format( i + 16, filepath, pot.get_izz(comm), \
                                 pot.s, pot.p, pot.d, pot.pnl, pot.dnl ))

  def _run_escan(self, comm, structure, norun):
    """ Runs escan only """
    from os.path import basename
    from ..opt import redirect

    self._write_incar(comm, structure, norun)
    if norun == False:
      with redirect(fout=self._cout(comm), ferr=self._cerr(comm), append=True) as oestreams: 
        assert comm.real, RuntimeError('Cannot run escan without mpi.')
        from ._escan import _call_escan
        _call_escan(comm)

  def _get_kpoint(self, structure, comm, norun):
    """ Returns deformed or undeformed kpoint. """
    from numpy import abs, sum, zeros, array, any, dot
    from numpy.linalg import inv
    from quantities import angstrom
    from ..physics import a0
    from ..crystal import to_voronoi
    # first get relaxed cell
    if norun == True: relaxed = structure.cell.copy()
    else:
      relaxed = zeros((3,3), dtype="float64")
      if comm.is_root:
        with open(self._POSCAR, "r") as file:
          file.readline() # number of atoms.
          # lattice vector by lattice vector
          for i in range(3): 
            relaxed[:,i] = array([float(u) for u in file.readline().split()[:3]])
        relaxed = relaxed / structure.scale * float(a0.rescale(angstrom))
      relaxed = comm.broadcast(relaxed)

    input = structure.cell 
    if (norun != True) and self.do_relax_kpoint and any(abs(relaxed-input) > 1e-12):
      kpoint = dot(dot(relaxed.T, inv(input.T)), self.kpoint)
    else:
      kpoint = self.kpoint
    return 1, kpoint[0], kpoint[1], kpoint[2], structure.scale / float(a0.rescale(angstrom))

  def copy(self, **kwargs):
    """ Performs deep-copy of object. 

	:Param kwargs: keyword arguments will overide the corresponding
           attribute. The value of the keyword argument is *not* deepcopied. 
           The attribute must exist. Attributes cannot be created here. 
    """
    from copy import deepcopy
    result = deepcopy(self)
    for key, value in kwargs.iteritems():
      assert hasattr(result, key),\
             ValueError("Attribute {0} does not exist and will not be created in copy.".format(key))
      setattr(result, key, value)
    return result

  def __setstate__(self, value):
    """ Takes care of different versions.
 
        Earlier versions may not have all the attributes that now exist by
        defaults. This routines adds them when unpickling.
    """
    self.__dict__.update(value)
    self.symlink = True
    for key, value in self.__class__().__dict__.iteritems():
       if not hasattr(self, key): setattr(self, key, value)
