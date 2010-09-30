""" A class for single-shot VASP calculation """
__docformat__ = "restructuredtext en"
import incar
from incar import Incar
from kpoints import Density
from ..opt.decorators import add_setter


class Launch(Incar):
  """ A class to launch a single vasp calculation """

  # program = "vasp.mpi"
  # """ shell command to launch vasp. """

  def __init__(self, inplace = True, workdir = None, species = None,
               kpoints = Density(), **kwargs ):
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
    from os import getcwd
    from ..opt import RelativeDirectory
    super(Launch, self).__init__() 

    self._workdir = RelativeDirectory(workdir)
    """ Filesystem location where temporary calculations are performed.  """
    self.species = species if species != None else {}
    """ Species in the system. """
    self.kpoints = kpoints
    """ kpoints for which to perform calculations. """
    self.inplace = inplace
    """ If True calculations are performed in the output directory. """
    self.print_from_all = False
    """ If True, will print from all nodes rather than just root. """

    # checks inplace vs workdir
    if self.inplace: 
      assert workdir == None, ValueError("Cannot use both workdir and inplace attributes.")

    # sets all other keywords as attributes.
    for key in kwargs.keys(): setattr(self, key, kwargs[key])

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

        - Writes INCAR file.
        - Writes KPOINTS file.
        - Writes POSCAR file.
        - Creates POTCAR file
        - Saves pickle of self.
    """
    import cPickle
    from copy import deepcopy
    from os.path import join, abspath
    from os import getcwd
    from shutil import copy
    from . import files, is_vasp_5
    from ..crystal import write_poscar, specie_list
    from ..opt.changedir import Changedir

    # creates incar file. Changedir makes sure that any calculations done to
    # obtain incar will happen in the tempdir. Only head node actually writes.
    # Also, to make sure in-place calculations will not end up with a bunch of
    # intermediate file, we copy self and make it non-inplace.
    this = deepcopy(self) if self.inplace else self
    this.inplace = False
    with Changedir(self._tempdir) as tmpdir:
      incar_lines = this.incar_lines(comm = comm)

    if comm.rank != 0: return # don't have any more business here.

    with open(join(self._tempdir, files.INCAR), "w") as incar_file: 
      incar_file.writelines(incar_lines)
  
    # creates kpoints file
    with open(join(self._tempdir, files.KPOINTS), "w") as kp_file: 
      kp_file.write( self.kpoints(self) if hasattr(self.kpoints, "__call__") \
                     else self.kpoints + "\n" )
  
    # creates poscar file
    with open(join(self._tempdir, files.POSCAR), "w") as poscar: 
      write_poscar(self._system, poscar, is_vasp_5)
  
    # creates POTCAR file
    with open(join(self._tempdir, files.POTCAR), 'w') as potcar:
      for s in specie_list(self._system):
        potcar.writelines( self.species[s].read_potcar() )

    path = join(abspath(self._tempdir), files.FUNCCAR)
    with Changedir(outdir) as outdir: # allows relative paths.
      with open(path, 'w') as file: cPickle.dump((self), file)

  def _run(self, comm):
     """ Isolates calls to vasp itself """
     from os.path import join
     from . import files
     from . import call_vasp
     from ..opt import redirect
     from ..opt.changedir import Changedir
     # moves to working dir only now.
     stdout = join(self._tempdir, files.STDOUT) 
     stderr = join(self._tempdir, files.STDERR) 
     is_notroot = False if comm == None else comm.rank != 0
     if is_notroot and self.print_from_all:
       stdout += ".{0}".format(comm.rank)
       stderr += ".{0}".format(comm.rank)
     elif is_notroot: 
       stdout = "/dev/null"
       stderr = "/dev/null"
     with Changedir(self._tempdir):
       with redirect(fout=stdout, ferr=stderr) as streams:
         if comm != None: comm.barrier()
         call_vasp(comm)
             
#    with open(join(self._tempdir, files.STDOUT), "w") as stdout:
#      with open(join(self._tempdir, files.STDERR), "w") as stderr:
#        vasp_proc = Popen( self.program, cwd = self._tempdir, \
#                           stdout = stdout, stderr = stderr, shell = True )
#        vasp_proc.wait() # Wait for completion. Could set a timer here.

  def _postrun(self, repat, outdir, comm, norun):
     """ Copies files back to outdir """
     from os.path import exists, join, isdir, realpath
     from os import makedirs
     from shutil import copy
     from . import files

     is_root = True if comm == None else comm.rank == 0

     if is_root: 
       if not exists(outdir): makedirs(outdir)
       if not exists(outdir): raise IOError, "%s directory does not exist." % (outdir)
       if not isdir(outdir):  raise IOError, "%s is not a directory." % (outdir)
     
     if comm != None: comm.barrier()


     notfound = []
     if realpath(self._tempdir) != realpath(outdir):
       for filename in set(repat).union(files.minimal):
         if not is_root: filename = str("%s.%i" % (filename, comm.rank))
         if exists( join(self._tempdir, filename) ): copy( join(self._tempdir, filename), outdir )
         elif is_root: notfound.append(filename)

     if not norun: assert len(notfound) == 0, IOError("Files %s were not found.\n" % (notfound))

  def __call__( self, structure=None, outdir = None, comm = None, repat = [], \
                norun = False, keep_tempdir=False):
    from os import getcwd
    from shutil import copy2 as copy
    from boost.mpi import world
    from ..opt import Tempdir, Changedir, RelativeDirectory

    # set up
    if structure != None: self._system = structure
    elif not hasattr(self, "_system"): raise RuntimeError, "Internal bug.\n"
    outdir = getcwd() if outdir == None else RelativeDirectory(outdir).path
    if comm == None: comm = world

    is_root = comm.rank == 0

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
      if not norun: self._run(comm)

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
    

       

