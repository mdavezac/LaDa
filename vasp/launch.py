""" A class for single-shot VASP calculation """
from incar import Incar
from kpoints import Density


class Launch(Incar):
  """ A class to launch a single vasp calculation """

  # program = "vasp.mpi"
  # """ shell command to launch vasp. """

  def __init__(self, workdir = None, indir = None, species = None, kpoints = Density(), **kwargs ):
    """ Initializes Launch instance.

        The initializer can take any number of keyword arguments. Those other
        than kpoints and species will be set as attributes: C{launch =
        Launch(whatnot = something)} is equivalent to
        >> launch = Launch()
        >> launch.whatnot = something
        This behavior is useful to pass non-standard incar arguments directly
        in initialization.
        @param workdir: working directory. Defaults to current directory at
           time of calculation.
        @param indir: input directory.
        @param species: Species in the system. 
        @type species: list of L{Specie}
        @param kpoints: Kpoint behavior.
        @type kpoints: see L{kpoints.Density} and L{kpoints.Gamma}
    """
    Incar.__init__(self) 

    self.workdir = workdir
    self.indir = indir
    # sets species
    if species != None: self.species = species
    self.kpoints =  kpoints

    # sets all other keywords as attributes.
    for key in kwargs.keys(): getattr(self, key).value = kwargs[key]

  def _find_species( self, structure ):
    """ Returns a list of species in the structure. 

        @raise AssertionError: when appropriate POTCAR or POTCAR.Z cannot be found.
    """
    from os.path import exists, join

    results = []
    for atom in structure.atoms:
      for s in self.species:
        if s.symbol == atom.type and not(s in results):
          results.append( s )
          if not exists(join( s.path, "POTCAR" )):
            raise AssertionError, "Could not find potcar in " + s.path
    return results

  def _prerun(self, comm):
    from os.path import join, exists
    from os import getcwd
    from shutil import copy
    from subprocess import Popen, PIPE
    from . import files
    from ..crystal import print_poscar
    from ..opt.changedir import Changedir

    # creates incar file
    with Changedir(self._tempdir) as tmpdir:
      incar_string = "".join((param.incar_string(self, comm) + "\n") for param in self)
    if comm != None:
      if comm.rank != 0: return # don't have any more business here.
    with open(join(self._tempdir, files.INCAR), "w") as incar_file: 
      incar_file.write(incar_string)
  
    # creates kpoints file
    with open(join(self._tempdir, files.KPOINTS), "w") as kpoints: 
      kpoints.write( self.kpoints(self) )
  
    # creates poscar file
    print_poscar(self._system, tuple(u.symbol for u in self.species), self._tempdir)
  
    # creates POTCAR file
    with open(join(self._tempdir, files.POTCAR), 'w') as potcar:
      for s in self._find_species(self._system):
        if not exists(join(s.path, "POTCAR")):  
          raise AssertionError, "Could not find potcar in " + s.path
        with open(join(s.path, "POTCAR"), "r") as infile: potcar.writelines(infile)

    # checks for CHGCAR
    indir = self.indir
    if indir == None: indir = getcwd()
    chgcar = join(indir, files.CHGCAR)
    if exists(chgcar): copy(chgcar, self._tempdir)

    # checks for WAVECAR
    wavecar = join(indir, files.WAVECAR) 
    if exists(chgcar): copy(chgcar, self._tempdir)

    # checks for EIGENVALUE
    eigenvalue = join(indir, files.EIGENVALUES) 
    if exists(eigenvalue): copy(eigenvalue, self._tempdir)

  def _run(self, comm):
     """ Isolates calls to vasp itself """
     from os.path import join
     from subprocess import Popen
     from . import files
     from . import call_vasp
     from ..opt import Redirect
     from ..opt.changedir import Changedir
     # moves to working dir only now.
     stdout = join(self._tempdir, files.STDOUT) 
     stderr = join(self._tempdir, files.STDERR) 
     if comm.rank != None:
       if comm.rank != 0:
         stdout += ".%i" % (comm.rank)
         stderr += ".%i" % (comm.rank)
     with Changedir(self._tempdir):
       with Redirect(Redirect.fortran.output, stdout) as stdout:
         with Redirect(Redirect.fortran.error, stderr) as stderr:
           if comm != None: comm.barrier()
           call_vasp(comm)
             
#    with open(join(self._tempdir, files.STDOUT), "w") as stdout:
#      with open(join(self._tempdir, files.STDERR), "w") as stderr:
#        vasp_proc = Popen( self.program, cwd = self._tempdir, \
#                           stdout = stdout, stderr = stderr, shell = True )
#        vasp_proc.wait() # Wait for completion. Could set a timer here.

  def _postrun(self, repat, outdir,  comm):
     """ Copies files back to outdir """
     from os.path import exists, join, isdir
     from os import makedirs
     from shutil import copy
     from . import files

     is_root = comm == None
     if comm != None: is_root = comm.rank == 0

     if is_root: 
       if not exists(outdir): makedirs(outdir)
       if not exists(outdir): raise IOError, "%s directory does not exist." % (outdir)
       if not isdir(outdir):  raise IOError, "%s is not a directory." % (outdir)
     
     if comm != None: comm.barrier()


     notfound = []
     for filename in set(repat).union(files.minimal):
       if not is_root: filename = str("%s.%i" % (filename, comm.rank))
       if exists( join(self._tempdir, filename) ): copy( join(self._tempdir, filename), outdir )
       elif is_root: notfound.append(filename)

     if len(notfound) != 0: raise IOError, "Files %s were not found.\n" % (notfound)

  def __call__(self, structure=None, outdir = None, comm = None, repat = []):
    from os.path import exists, join, abspath, expanduser
    from os import getcwd
    from shutil import copy2 as copy
    from ..opt.tempdir import Tempdir

    # set up
    if structure != None: self._system = structure
    elif not hasattr(self, "_system"): raise RuntimeError, "Internal bug.\n"
    if outdir == None: outdir = self.indir
    if outdir == None: outdir = getcwd()
    outdir = abspath(expanduser(outdir))

    is_root = comm == None
    if comm != None: is_root = comm.rank == 0 

    # creates temporary working directory
    workdir = self.workdir
    if workdir == None: workdir = getcwd()
    with Tempdir(workdir=workdir, comm=comm) as self._tempdir: 
      # We do not move to working directory to make copying of files from indir
      # or outdir (as relative paths) possible.
      # creates INCAR and KPOINTS.
      # copies POSCAR, POTCAR, possibly WAVECAR and CHGCAR from indir
      self._prerun(comm)

      # runs vasp.
      self._run(comm)

      # now copies data back
      self._postrun(repat, outdir, comm)

    # deletes system attribute.
    del self._system

