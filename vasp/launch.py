""" A class for single-shot VASP calculation """
__docformat__ = "restructuredtext en"
__all__ = ['Launch']
from incar import Incar
from kpoints import Density
from ..opt.decorators import add_setter


class Launch(Incar):
  """ A class to launch a single vasp calculation """

  def __init__( self, inplace = True, workdir = None, species = None, \
                kpoints = None, **kwargs ):
    """ Initializes Launch instance.

        :Parameters:
          workdir : str
            working directory. Defaults to output directory.
          inplace : boolean
            In-place calculations. Eg no workdir.
          species
           Species in the system. It is a dictionary of `vasp.Specie`.
           Individual items can be set using `Launch.add_specie`.
          kpoints : str or callable
            A string describing the kpoint mesh (in VASP's KPOINT format), or a
            callable returning such a string.
    """
    from ..opt import RelativeDirectory
    from .. import launch_vasp_as_library, vasp_library, vasp_program
    super(Launch, self).__init__() 

    self._workdir = RelativeDirectory(workdir)
    """ Filesystem location where temporary calculations are performed.  """
    self.species = species if species is not None else {}
    """ Species in the system. """
    self.kpoints = '\n0\nM\n2 2 2\n0 0 0\n' if kpoints is None else kpoints 
    """ kpoints for which to perform calculations. """
    self.inplace = inplace
    """ If True calculations are performed in the output directory. """
    self.print_from_all = False
    """ If True, will print from all nodes rather than just root. """
    self.symlink = kwargs.pop('symlink', False)
    """ If True, prefers symlinks where possible. """
    self.vasp_library = vasp_library
    """ Path to vasp-library. 
    
        If None, then use global default `lada.vasp_library`.
    """
    self.program = vasp_program
    """ Name/fullpath of vasp program. """
    self.launch_as_library = launch_vasp_as_library
    """ Whether to launch vasp as library(True) or a program(False). """

    # checks inplace vs workdir
    if self.inplace: 
      assert workdir is None, ValueError("Cannot use both workdir and inplace attributes.")

  @property
  def workdir(self):
    """ Directory where calculations are performed. 

        If `Launch.inplace` is true, then this parameter is ignored, and
        calculations are performed in the output directory.
    """
    return self._workdir.path
  @workdir.setter
  def workdir(self, value): self._workdir.path = value

  def _prerun(self, comm, outdir):
    """ Sets things up prior to calling VASP. 

        Performs the following actions.

        - Writes POSCAR file.
        - Writes INCAR file.
        - Writes KPOINTS file.
        - Creates POTCAR file
        - Saves pickle of self.
    """
    import cPickle
    from copy import deepcopy
    from os.path import join, abspath
    from . import files, is_vasp_5
    from ..crystal import write_poscar, specie_list
    from ..opt.changedir import Changedir

    # creates poscar file. Might be overwriten by restart.
    if comm.is_root:
      with open(join(self._tempdir, files.POSCAR), "w") as poscar: 
        which_type = is_vasp_5(self.vasp_library) if self.launch_as_library else False
        write_poscar(self._system, poscar, which_type)

    # creates incar file. Changedir makes sure that any calculations done to
    # obtain incar will happen in the tempdir. Only head node actually writes.
    # Also, to make sure in-place calculations will not end up with a bunch of
    # intermediate file, we copy self and make it non-inplace.
    this = deepcopy(self) if self.inplace else self
    this.inplace = False
    with Changedir(self._tempdir) as tmpdir:
      incar_lines = this.incar_lines(comm = comm)

    if not comm.is_root: return # don't have any more business here.


    # creates INCAR file. Note that POSCAR file might be overwritten here by Restart.
    with open(join(self._tempdir, files.INCAR), "w") as incar_file: 
      incar_file.writelines(incar_lines)
  
    # creates kpoints file
    with open(join(self._tempdir, files.KPOINTS), "w") as kp_file: 
      self.write_kpoints(kp_file)
  
    # creates POTCAR file
    with open(join(self._tempdir, files.POTCAR), 'w') as potcar:
      for s in specie_list(self._system):
        potcar.writelines( self.species[s].read_potcar() )

    path = join(abspath(self._tempdir), files.FUNCCAR)
    with Changedir(outdir) as outdir: # allows relative paths.
      with open(path, 'w') as file: cPickle.dump(self, file)

  def _run(self, comm, minversion):
     """ Isolates calls to vasp itself """
     from os.path import join
     from . import files
     from . import call_vasp
     from ..opt import redirect, which
     from ..opt.changedir import Changedir
     # moves to working dir only now.
     stdout = join(self._tempdir, files.STDOUT) 
     stderr = join(self._tempdir, files.STDERR) 
     if (not comm.is_root) and getattr(self, "print_from_all", False):
       stdout += ".{0}".format(comm.rank)
       stderr += ".{0}".format(comm.rank)
     elif not comm.is_root: 
       stdout = "/dev/null"
       stderr = "/dev/null"
     with Changedir(self._tempdir):
       if self.launch_as_library: 
         with redirect(fout=stdout, ferr=stderr) as streams:
           assert comm.real, ValueError("Cannot call vasp without mpi.")
           call_vasp(self.vasp_library, comm, minversion)
       else:
         try: program = which(self.program)
         except: program = self.program
         comm.external(program, out=stdout, err=stderr)
             
  def _postrun(self, repat, outdir, comm, norun):
     """ Copies files back to outdir """
     from os.path import exists, join, isdir, realpath
     from os import makedirs
     from shutil import copy
     from . import files

     if comm.is_root: 
       if not exists(outdir): makedirs(outdir)
       if not exists(outdir): raise IOError, "%s directory does not exist." % (outdir)
       if not isdir(outdir):  raise IOError, "%s is not a directory." % (outdir)
       # Appends INCAR and CONTCAR to OUTCAR:
       with open(files.OUTCAR, 'a') as outcar:
         outcar.write('\n################ INCAR ################\n')
         with open(files.INCAR, 'r') as incar: outcar.write(incar.read())
         outcar.write('\n################ END INCAR ################\n')
         outcar.write('\n################ CONTCAR ################\n')
         if self.isif == 0 or self.nsw == 0 or self.ibrion == -1:
           with open(files.POSCAR, 'r') as contcar: outcar.write(contcar.read())
         else:
           with open(files.CONTCAR, 'r') as contcar: outcar.write(contcar.read())
         outcar.write('\n################ END CONTCAR ################\n')

     
     comm.barrier()


     notfound = []
     if realpath(self._tempdir) != realpath(outdir):
       for filename in set(repat).union(files.minimal):
         if not comm.is_root: filename = str("%s.%i" % (filename, comm.rank))
         if exists( join(self._tempdir, filename) ): copy( join(self._tempdir, filename), outdir )
         elif comm.is_root: notfound.append(filename)

     if not norun: assert len(notfound) == 0, IOError("Files %s were not found.\n" % (notfound))

  def __call__( self, structure=None, outdir = None, comm = None, repat = None, \
                norun = False, keep_tempdir=False, minversion=0):
    from os import getcwd
    from ..opt import Tempdir, Changedir, RelativeDirectory
    from ..mpi import Communicator

    # set up
    if structure is not None: self._system = structure
    elif not hasattr(self, "_system"): raise RuntimeError, "Internal bug.\n"
    outdir = getcwd() if outdir is None else RelativeDirectory(outdir).path
    comm = Communicator(comm, with_world=True)
    repat = [] if repat is None else repat

    # creates temporary working directory
    if self.inplace: context = Changedir(outdir, comm=comm) 
    else:            context = Tempdir(workdir = self.workdir, comm = comm)
    with context as self._tempdir:
      # We do not move to working directory to make copying of files from indir
      # or outdir (as relative paths) possible.
      # creates INCAR and KPOINTS.
      # copies POSCAR, POTCAR, possibly WAVECAR and CHGCAR from indir
      self._prerun(comm, outdir)

      # runs vasp.
      if not norun: self._run(comm, minversion)

      # now copies data back
      self._postrun(repat, outdir, comm, norun)

    # deletes system attribute.
    del self._system


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
    Incar.__setstate__(self, args)
    for key, value in self.__class__().__dict__.iteritems():
       if not hasattr(self, key): setattr(self, key, value)

