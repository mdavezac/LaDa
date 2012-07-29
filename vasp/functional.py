""" Sub-package containing the functional. """
__docformat__ = "restructuredtext en"
__all__ = ['Vasp']
from ..functools import stateless, assign_attributes
from ..functools.block import AttrBlock
from ..misc import add_setter
from extract import Extract


class Vasp(AttrBlock):
  """ Interface to VASP code. """
  Extract = staticmethod(Extract)
  """ Extraction class. """

  def __init__(self, copy=None, species=None, kpoints=None, **kwargs):
    """ Initializes vasp class. """
    from .keywords import BoolKeyword, Magmom, System, Npar, ExtraElectron,    \
                          NElect, Algo, Ediff, Ediffg, Encut, EncutGW, IStart, \
                          ICharge, IStruc, UParams, PrecFock, Precision, Nsw,  \
                          Isif, Ibrion, Relaxation
    from ..functools import TypedKeyword, from ..functools import ChoiceKeyword
    super(Vasp, self).__init__()

    # copies values from other functional.
    if copy is not None: 
      self.params.update(copy.params)
      self.special.update(copy.special)
      for key, value in copy.__dict__.iteritems():
        if key in kwargs: continue
        elif key == 'params': continue
        elif key == 'special': continue
        elif hasattr(self, key): setattr(self, key, value)

    self.species = species if species is not None else {}
    """ Species in the system. """
    self.kpoints = kpoints if kpoints is not None \
                   else "\n0\nM\n4 4 4\n0 0 0"
    """ kpoints for which to perform calculations. """

    self.program = kwargs.get('program', None)
    """ Path to vasp program. 
    
        Can be one of the following:

          - None: defaults to :py:attr:`~lada.vasp_program`.
            :py:attr:`~lada.vasp_program` can take the same values as described
            here, except for None.
          - string: Should be the path to the vasp executable. It can be either
            a full path, or an executable within the envirnoment's $PATH
            variable.
          - callable: The callable is called with a :py:class:`~lada.vasp.Vasp`
            as sole argument. It should return a string, as described above.
            In other words, different vasp executables can be used depending on
            the parameters. 
    """
    self.addgrid = BoolKeyword()
    """ Adds additional support grid for augmentation charge evaluation. 

        Can be only True or False (or None for VASP_ default).

        .. seealso:: 

           ADDGRID_

           .. _ADDGRID: http://cms.mpi.univie.ac.at/wiki/index.php/ADDGRID
    """
    self.ispin   = ChoiceKeyword(values=(1, 2))
    """ Whether to perform spin-polarized or spin-unpolarized calculations.

        Can be only 1 or 2 (or None for VASP_ default).

        .. seealso:: 

           ISPIN_ 

           .. _ISPIN: http://cms.mpi.univie.ac.at/wiki/index.php/ISPIN
    """
    self.istart    = ChoiceKeyword(values=range(5))
    """ Starting point of calculation.

        Can take a value between 1 and 4 included (and None for VASP_ default).

        .. seealso:: 

           ISTART_ 

           .. _ISTART: http://cms.mpi.univie.ac.at/wiki/index.php/ISTART
    """
    self.isym      = ChoiceKeyword(values=range(3))
    """ Symmetry scheme.

        .. seealso:: 
         
           ISYM_

           .. _ISYM: http://cms.mpi.univie.ac.at/vasp/guide/node115.html
    """ 
    self.lmaxmix   = TypedKeyword(type=int)
    """ Cutoff *l*-quantum number of PAW charge densities passed to mixer 

        .. seealso:: 

           LMAXMIX_ 

           .. _LMAXMIX: http://cms.mpi.univie.ac.at/wiki/index.php/LMAXMIX
    """
    self.lorbit    = ChoiceKeyword(values=(0, 1, 2, 5, 10, 11, 12))
    """ Decides whether PROOUT and PROOCAR are writtent to disk.

        Can be one of 0|1|2|5|10|11|12|None. 

        .. seealso:: 

           LORBIT_ 

           .. _LORBIT: http://cms.mpi.univie.ac.at/wiki/index.php/LORBIT
    """
    self.nbands    = TypedKeyword(type=int)
    self.nomega    = TypedKeyword(type=int)
    self.nupdown   = TypedKeyword(type=int)
    self.symprec   = TypedKeyword(type=float)
    self.lwave     = BoolKeyword(value=False)
    self.lcharg    = BoolKeyword(value=True)
    self.lvtot     = BoolKeyword(value=False)
    self.lrpa      = BoolKeyword()
    self.loptics   = BoolKeyword()
    self.lpead     = BoolKeyword()
    self.nelm      = TypedKeyword(type=int)
    self.nelmin    = TypedKeyword(type=int)
    self.nelmdl    = TypedKeyword(type=int)
    self.ngx       = TypedKeyword(type=int)
    self.ngy       = TypedKeyword(type=int)
    self.ngz       = TypedKeyword(type=int)
    self.nonscf    = BoolKeyword()
    """ If True, performs a non-self consistent calculation.

        The value of this keyword is checked by :py:attr:`icharg` and used
        appropriately.
    """
    
    self.magmom    = Magmom()
    """ Sets the initial magnetic moments on each atom.
    
        There are three types of usage: 
    
        - if None or False, does nothing
        - if calculations are not spin-polarized, does nothing.
        - if a string, uses that as for the MAGMOM_ keyword
        - if True and at least one atom in the structure has a non-zero
          ``magmom`` attribute, then creates the relevant moment input for VASP_
    
        If the calculation is **not** spin-polarized, then the magnetic moment
        tag is not set.
    
        .. note:: Please set by hand for non-collinear calculations

        .. seealso:: MAGMOM_

        .. _MAGMOM: http://cms.mpi.univie.ac.at/wiki/index.php/MAGMOM
    """
    self.system    = System()
    """ System title to use for calculation.
    
        - If None and ... 
           - if the structure has a ``name`` attribute, uses that as the
             calculations title
           - else does not use SYSTEM_ tag
        - If something else which is convertable to a string,  and ...
           - if the structure has a ``name`` attribute, uses 'string: name' as
             the title
           - otherwise, uses the string
    
        .. seealso:: SYSTEM_
    
        .. _SYSTEM: http://cms.mpi.univie.ac.at/vasp/guide/node94.html>
    """
    self.npar      = Npar()
    """ Parallelization over bands. 
    
        Npar defines how many nodes work on one band.
        It can be set to a particular number:
    
        >>> vasp.npar = 2
    
        Or it can be deduced automatically. Different schemes are available:
        
          - power of two: npar is set to the largest power of 2 which divides the
            number of processors.
   
            >>> vasp.npar = "power of two"
    
            If the number of processors is not a power of two, prints nothing.
    
          - square root: npar is set to the square root of the number of processors.
    
            >>> vasp.npar = "sqrt"
        
    
        .. seealso: `NPAR <http://cms.mpi.univie.ac.at/vasp/guide/node138.html>`_
    """
    self.extraelectron = ExtraElectron()
    """ Number of electrons relative to neutral system.
        
        Gets the number of electrons in the (neutral) system. Then adds value to
        it and computes with the resulting number of electrons.
    
        >>> vasp.extraelectron =  0  # charge neutral system
        >>> vasp.extraelectron =  1  # charge -1 (1 extra electron)
        >>> vasp.extraelectron = -1  # charge +1 (1 extra hole)
    
        .. seealso:: `NELECT <http://cms.mpi.univie.ac.at/wiki/index.php/NELECT>`_
    """
    self.nelect = NElect()
    """ Sets the absolute number of electrons.
        
        Disables :py:attr:`lada.vasp.functional.Functional.extraelectron` if set to
        something other than None.
    
        .. seealso:: `NELECT <http://cms.mpi.univie.ac.at/wiki/index.php/NELECT>`_
    """
    self.algo = Algo()
    """ Electronic minimization. 
    
        Defines the kind of algorithm vasp will run.
          - very fast
          - fast, f (default)
          - normal, n
          - all, a
          - damped, d 
          - Diag 
          - conjugate, c (vasp 5)
          - subrot (vasp 5)
          - eigenval (vasp 5)
          - Nothing (vasp 5)
          - Exact  (vasp 5)
          - chi
          - gw
          - gw0
          - scgw
          - scgw0
    
        If :py:data:`is_vasp_4 <lada.is_vasp_4>` is an existing configuration
        variable of :py:mod:`lada` the parameters marked as vasp 5 will fail.
    
        .. warning:: The string None is not  allowed, as it would lead to
           confusion with the python object None. Please use "Nothing" instead.
           The python object None will simply not print the ALGO keyword to the
           INCAR file.
    
        .. note:: By special request, "fast" is the default algorithm.
    
        .. seealso:: `ALGO <http://cms.mpi.univie.ac.at/vasp/vasp/ALGO_tag.html>`_
    """ 
    self.ediff = Ediff()
    """ Sets the convergence criteria (per atom) for electronic minimization.
    
        - value > 0e0: the tolerance is multiplied by the number of atoms in the
          system. This makes tolerance consistent from one system to the next.
        - value < 0e0: tolerance is given as absolute value, without multiplying
          by size of system.
    
        .. seealso:: `EDIFF <http://cms.mpi.univie.ac.at/vasp/guide/node105.html>`_
    """
    self.ediffg = Ediffg()
    """ Sets the convergence criteria (per atom) for ionic minimization.
    
        - value > 0e0: the tolerance is multiplied by the number of atoms in the
          system. This makes tolerance consistent from one system to the next.
        - value < 0e0: tolerance is given as is (negative), and applies to forces.
    
        .. seealso:: `EDIFFG <http://cms.mpi.univie.ac.at/vasp/guide/node107.html>`_
    """
    self.encut = Encut()
    """ Defines cutoff factor for calculation. 
    
        There are three ways to set this parameter:
    
        - if value is floating point and 0 < value <= 3: then the cutoff is
          ``value * ENMAX``, where ENMAX is the maximum recommended cutoff for
          the species in the system.
        - if value > 3 eV, then prints encut is exactly value (in eV). Any energy
          unit is acceptable.
        - if value < 0 eV or None, does not print anything to INCAR. 
        
        .. seealso:: `ENCUT <http://cms.mpi.univie.ac.at/vasp/vasp/ENCUT_tag.html>`_
    """
    self.encutgw = EncutGW()
    """ Defines cutoff factor for GW calculation. 
  
        There are three ways to set this parameter:
  
        - if value is floating point and 0 < value <= 3: then the cutoff is
          ``value * ENMAX``, where ENMAX is the maximum recommended cutoff for
          the species in the system.
        - if value > 3 eV, then prints encut is exactly value (in eV). Any energy
          unit is acceptable.
        - if value < 0 eV or None, does not print anything to INCAR. 
        
        .. seealso:: `ENCUTGW
          <http://cms.mpi.univie.ac.at/vasp/vasp/ENCUTGW_energy_cutoff_response_function.html>`_
    """
    self.istart = IStart()
    """ Starting wavefunctions.
    
        It is best to keep this attribute set to -1, in which case, LaDa takes
        care of copying the relevant files.
    
          - -1: Automatically determined by LaDA. Depends on the value of
                restart_ and the existence of the relevant files.
    
          - 0: Start from scratch.
    
          - 1: Restart with constant cutoff.
    
          - 2: Restart with constant basis.
    
          - 3: Full restart, including TMPCAR.
    
        .. note::
        
           Files are copied right before the calculation takes place, not
           before.
    
        .. seealso:: ISTART_
    
        .. _ISTART: http://cms.mpi.univie.ac.at/wiki/index.php/ISTART
        .. _restart: :py:attr:`~lada.vasp.functional.Functional.restart`
    """ 
    self.icharg = ICharge()
    """ Charge from which to start. 
    
        It is best to keep this attribute set to -1, in which case, LaDa takes
        care of copying the relevant files.
    
          - -1: Automatically determined by LaDA. Depends on the value of
                restart_ and the existence of the relevant files. Also takes
                care of non-scf bit.
    
          - 0: Tries to restart from wavefunctions. Uses the latest WAVECAR
               file between the one currently in the output directory and the
               one in the restart directory (if speciefied). Sets nonscf_ to
               False.
    
               .. note:: CHGCAR is also copied, just in case.
    
          - 1: Tries to restart from wavefunctions. Uses the latest WAVECAR
               file between the one currently in the output directory and the
               one in the restart directory (if speciefied). Sets nonscf_ to
               False.
    
          - 2: Superimposition of atomic charge densities. Sets nonscf_ to
               False.
    
          - 4: Reads potential from POT file (VASP-5.1 only). The POT file is
               deduced the same way as for CHGAR and WAVECAR above.  Sets
               nonscf_ to False.
    
          - 10, 11, 12: Same as 0, 1, 2 above, but also sets nonscf_ to True.
               This is a shortcut. The value is actually kept to 0, 1, or 2:
    
               >>> vasp.icharg = 10
               >>> vasp.nonscf, vasp.icharg
               (True, 0)
    
        .. note::
        
           Files are copied right before the calculation takes place, not before.
    
        .. seealso:: ICHARG_
    
        .. _ICHARG: http://cms.mpi.univie.ac.at/wiki/index.php/ICHARG
        .. _restart: :py:attr:`~lada.vasp.functional.Functional.restart`
        .. _nonscf: :py:attr:`~lada.vasp.functional.Functional.nonscf`
    """ 
    self.functional.istruc = IStruc()
    """ Initial structure. 
    
        Determines which structure is written to the POSCAR. In practice, it
        makes it possible to restart a crashed job from the latest contcar.
        There are two possible options:
  
          - auto: LaDa determines automatically what to use. If a CONTCAR
                  exists in either the current directory or in the restart
                  directory (if any), then uses the latest. Otherwise, uses
                  input structure.
          - scratch: Always uses input structure.
  
        If the run was given the ``overwrite`` option, then always uses the
        input structure.

        .. note:: There is no VASP equivalent to this option.
    """
    self.ldauprint = UParams()
    """ Sets U, nlep, and enlep parameters. 
   
        The U, nlep, and enlep parameters of the atomic species are set at the
        same time as the pseudo-potentials. This object merely sets up the incar
        with right input.
    
        However, it does accept one parameter, which can be "off", "on", "occ" or
        "all", and defines the level of verbosity of VASP (with respect to U and nlep).
    
        .. seealso:: `LDAU, LDAUTYPE, LDAUL, LDAUPRINT
          <http://cms.mpi.univie.ac.at/vasp/vasp/On_site_Coulomb_interaction_L_S_DA_U.html>`_
    """
    self.precfock = PrecFock()
    """ Sets up FFT grid in hartree-fock related routines.
        
        Allowable options are:
    
        - low
        - medium
        - fast
        - normal
        - accurate
    
        .. seealso:: PRECFOCK_
        
        .. _PRECFOCK: http://cms.mpi.univie.ac.at/wiki/index.php/PRECFOCK
    """
    self.precision = Precision()
    """ Sets accuracy of calculation. 
    
        - accurate (default)
        - low
        - medium
        - high
        - single
    
        .. seealso:: PREC_
        
        .. _PREC: http://cms.mpi.univie.ac.at/wiki/index.php/PREC
    """
    self.nsw = Nsw()
    """ Maxium number of ionic iterations. 

        .. seealso:: NSW_

        .. _NSW: http://cms.mpi.univie.ac.at/wiki/index.php/NSW
    """
    self.ibrion = Ibrion()
    """ Ions/cell-shape/volume optimization method.
    
        Can only take a restricted set of values: -1 | 0 | 1 | 2 | 3 | 5 | 6 | 7 | 8 | 44.

        .. seealso:: IBRION_

        .. _IBRION: cms.mpi.univie.ac.at/wiki/index.php/IBRIONN
    """
    self.isif = Isif()
    """ Degree of librerty to optimize during geometry optimization

        .. seealso:: ISIF_

        .. _ISIF: http://cms.mpi.univie.ac.at/vasp/guide/node112.html
    """
    self.relaxation = Relaxation()
    """ Short-cut for setting up relaxation. 

        It accepts two parameters:
        
          - static: for calculation without geometric relaxation.
          - combination of ionic, volume, cellshape: for the type of relaxation
            requested.

        It makes sure that isif_, ibrion_, and nsw_ take the right value for the
        kind of relaxation.
    """




    # sets all known keywords as attributes.
    for key, value in kwargs.iteritems():
      if hasattr(self, key): setattr(self, key, value)

  def __call__( self, structure, outdir=None, comm=None, overwrite=False, 
                **kwargs):
    """ Calls vasp program. """
    result = None
    for program in self.iter(structure, outdir=outdir, comm=comm, overwrite=overwrite, **kwargs):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False):
        result = program
        continue
      # otherwise, it should yield a Program tuple to execute.
      program.start(comm)
      program.wait()
    # Last yield should be an extraction object.
    if not result.success:
      raise RuntimeError("Vasp failed to execute correctly.")
    return result

  @assign_attributes(ignore=['overwrite', 'comm'])
  @stateless
  def iter(self, structure, outdir=None, comm=None, overwrite=False, **kwargs):
    """ Performs a vasp calculation 
     
        If successfull results (see :py:attr:`extract.Extract.success`) already
        exist in outdir, calculations are not repeated. Instead, an extraction
        object for the stored results are given.

        :param structure:  
            :py:class:`~lada.crystal.Structure` structure to compute, *unless*
            a CONTCAR already exists in ``outdir``, in which case this
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
            If False, will check whether a successfull calculation exists. If
            one does, then does not execute. 
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
    from ..process.program import ProgramProcess

    # check for pre-existing and successfull run.
    if not overwrite:
      extract = self.Extract(outdir)
      if extract.success:
        yield extract # in which case, returns extraction object.
        return
    
    # makes functor stateless/reads structure from CONTCAR if requested and appropriate.
    if kwargs.pop("restart_from_contcar", self.restart_from_contcar): 
      path = join(outdir, CONTCAR)
      if exists(path):
        try: contstruct = read.poscar(path, list(specieset(structure)))
        except: pass
        else:
          # copies poscar info to structure.
          # this function should be stateless at this point, so it does not
          # matter that we change structure.
          for a, b in zip(structure, contstruct):
            a.pos, a.type = b.pos, b.type
          structure.cell = contstruct.cell
          structure.scale = contstruct.scale
    if len(structure) == 0: raise ValueError("Structure is empty.")
    if abs(det(structure.cell)) < 1e-8: raise ValueError("Structure with zero volume.")
    if abs(structure.scale) < 1e-8: raise ValueError("Structure with null scale.")

    # copies/creates file environment for calculation.
    self.bringup(structure, outdir, comm)
    # figures out what program to call.
    program = self.program if self.program is not None else vasp_program
    if hasattr(program, '__call__'): program = program(self)
    # creates a process with a callback to bring-down environment once it is
    # done.
    def onfinish(process, error):  self.bringdown(outdir, structure)
    yield ProgramProcess( program, cmdline=[], outdir=outdir,
                          onfinish=onfinish, stdout='stdout', stderr='stderr',
                          dompi=True )
    # yields final extraction object.
    yield Extract(outdir)

  def bringup(self, structure, outdir, comm):
    """ Creates all input files necessary to run results.

        Performs the following actions.

        - Writes POSCAR file.
        - Writes INCAR file.
        - Writes KPOINTS file.
        - Creates POTCAR file
        - Saves pickle of self.
    """
    import cPickle
    from ..crystal import specieset
    from ..misc.changedir import Changedir
    from . import files

    with Changedir(outdir) as tmpdir:
      # creates INCAR file. Note that POSCAR file might be overwritten here by Restart.
      self.write_incar(structure, comm=comm)
  
      # creates kpoints file
      with open(files.KPOINTS, "w") as kp_file: 
        self.write_kpoints(kp_file, structure)
  
      # creates POTCAR file
      with open(files.POTCAR, 'w') as potcar:
        for s in specieset(structure):
          potcar.writelines( self.species[s].read_potcar() )
    
      with open(files.FUNCCAR, 'w') as file:
        cPickle.dump(self, file)
    

  def bringdown(self, directory, structure):
     """ Copies contcar to outcar. """
     from . import files
     from ..misc import Changedir

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
         outcar.write("""from {0.__class__.__module__} import {0.__class__.__name__}\n"""\
                      """structure = {1}\n"""\
                      .format(structure, repr(structure).replace('\n', '\n            ')))
         outcar.write('\n################ END INITIAL STRUCTURE ################\n')
         outcar.write('\n################ FUNCTIONAL ################\n')
         outcar.write(repr(self))
         outcar.write('\n################ END FUNCTIONAL ################\n')


  def write_incar(self, structure, path=None, comm=None):
    """ Writes incar file. """
    from lada import default_comm
    from ..misc import RelativePath
    from .files import INCAR

    # check what type path is.
    # if not a file, opens one an does recurrent call.
    if path is None: path = INCAR
    if not hasattr(path, "write"):
      with open(RelativePath(path).path, "w") as file:
        self.write_incar(structure, path=file, comm=comm)
      return

    if comm is None: comm = default_comm
    map = self.input_map(structure=structure, vasp=self, comm=comm)
    length = max(len(u) for u in map)
    for key, value in map.iteritems():
      path.write('{0: >{length}} = {1}\n'.format(key, value, length=length))

  def write_kpoints(self, file, structure, kpoints=None):
    """ Writes kpoints to a stream. """
    if kpoints == None: kpoints = self.kpoints
    if isinstance(self.kpoints, str): file.write(self.kpoints)
    elif hasattr(self.kpoints, "__call__"):
      self.write_kpoints(file, structure, self.kpoints(self, structure))
    else: # numpy array or such.
      file.write("Explicit list of kpoints.\n{0}\nCartesian\n".format(len(self.kpoints)))
      for kpoint in self.kpoints:
        file.write("{0[0]} {0[1]} {0[2]} {1}\n".format(kpoint, 1 if len(kpoint) == 3 else kpoint[3]))

  def __repr__(self, skip=None):
    """ Returns a python script representing this object. """
    from .. import verbose_representation
    if skip is None: skip = not verbose_representation

    # creates a default vasp instance to compare to.
    compare = self.__class__()
    params = compare.params.keys()

    # will hold classes from modules.
    modules = {}
    modules[self.__class__.__module__] = [self.__class__.__name__]
    addparam = {}
    noaddparam = {}
    special = {}
    # now gather vasp parameters and check their length.
    for name, value in self.params.items():
      if skip and value is None: continue
      if name in params:
        if skip and value is None: continue
        try: # check if is default.
          if skip and value == compare.params[name]: continue
        except: pass
        noaddparam[name] = len(name), value
      else:
        if value is None: continue
        addparam[name] = len(name), value
        module = value.__class__.__module__ 
        classname = value.__class__.__name__ 
        if module in modules: modules[module].append(classname)
        else: modules[module] = [classname]
    # if a special parameter, then is non-default.
    for name, value in self.special.items():
      if skip and value.value is None: continue
      try: # check if is default.
        if skip and value.__class__ is compare.special[name].__class__ \
           and getattr(self, name) == getattr(compare, name): continue
      except: pass
      try: 
        if value.__class__ is compare.special[name].__class__:
          noaddparam[name] = len(name), value.value
          continue
      except: pass
      special[name] = len(name), value
      module = value.__class__.__module__ 
      classname = value.__class__.__name__ 
      if module in modules: modules[module].append(classname)
      else: modules[module] = [classname]
    # adds kpoints
    if hasattr(self.kpoints, "__call__"):
      # checks for user module.
      module = self.kpoints.__class__.__module__ 
      classname = self.kpoints.__class__.__name__ 
      if module in modules: modules[module].append(classname)
      else: modules[module] = [classname]
      noaddparam['kpoints'] = len('kpoints'), self.kpoints
    else:
      try: 
        if skip == False or self.kpoints != compare.kpoints:
          noaddparam['kpoints'] = len('kpoints'), self.kpoints
      except: noaddparam['kpoints'] = len('kpoints'), self.kpoints


    if not self.restart_from_contcar: 
      noaddparam['restart_from_contcar'] = len('restart_from_contcar'), False
    # adds objects in __dict__
    class Dummy: pass
    for key, value in self.__dict__.iteritems():
      if key == 'special' or key == 'params' or key == 'species': continue
      try: 
        if skip and getattr(compare, key, Dummy) == value: continue
      except: pass
      noaddparam[key] = len(key), value
      module = value.__class__.__module__
      classname = value.__class__.__name__ 
      if module == '__builtin__': continue
      if module in modules: modules[module].append(classname)
      else: modules[module] = [classname]

    # now write stuff to string.
    string = "functional = {0.__class__.__name__}()\n".format(self)

    def sortme(a): return a.lower()
    l = [k for k in noaddparam.iterkeys() if k[0] != '_']
    if len(l) != 0:
      length = max([len(k) for k in l])
      for key in sorted(l, key=sortme):
        string += 'functional.{0:<{length}} = {1[1]!r}\n'.format(key, noaddparam[key], length=length)
    l = [k for k in addparam.iterkeys() if k[0] != '_']
    if len(l) != 0:
      length = max([len(k) for k in l])
      for key in sorted(l, key=sortme):
        string += 'functional.add_param = {0!r: >{length}}, {1[1]!r}\n'\
                  .format(key, addparam[key], length=length)
    l = [k for k in noaddparam.iterkeys() if k[0] == '_']
    if len(l) != 0:
      length = max([len(k) for k in l])
      for key in sorted(l, key=sortme):
        string += 'functional.{0:<{length}} = {1[1]!r}\n'.format(key, noaddparam[key], length=length)
    l = [k for k in addparam.iterkeys() if k[0] == '_']
    if len(l) != 0:
      length = max([len(k) for k in l])
      for key in sorted(l, key=sortme):
        string += 'functional.add_param = {0!r: >{length}}, {1[1]!r}\n'\
                  .format(key, addparam[key], length=length)
    l = [k for k in special.iterkeys() if k[0] != '_']
    if len(l) != 0:
      length = max([len(k) for k in l])
      for key in sorted(l, key=sortme):
        string += 'functional.{0:<{length}} = {1[1]!r}\n'.format(key, special[key], length=length)
    l = [k for k in special.iterkeys() if k[0] == '_']
    if len(l) != 0:
      length = max([len(k) for k in l])
      for key in sorted(l, key=sortme):
        string += 'functional.add_param = {0!r: >{length}}, {1[1]!r}\n'\
                  .format(key, special[key], length=length)

    # adds species.
    if len(self.species) > 0:
      length = max([len(k) for k in self.species])
      for name, specie in self.species.items():
        string += "functional.species[{0!r}] {2:<{3}}= {1!r}\n".format(name, specie, '', length-len(name))
        module = specie.__class__.__module__ 
        classname = specie.__class__.__name__ 
        if module in modules: modules[module].append(classname)
        else: modules[module] = [classname]

    # adds user modules above repr string.
    header = ""
    for name in sorted(modules.keys()):
      mods = list(set(modules[name]))
      header += "from {0} import {1}".format(name, mods[0])
      for v in mods[1:]: header += ", {0}".format(v)
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
    super(Vasp, self).__setstate__(args)
    for key, value in self.__class__().__dict__.iteritems():
       if not hasattr(self, key): setattr(self, key, value)

