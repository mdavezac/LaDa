""" A class for single-shot VASP calculation """
from incar import Incar
from kpoints import Density


class Launch(Incar):
  """ A class to launch a single vasp calculation """

  program = "vasp.mpi"
  """ shell command to launch vasp. """

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
    for key in kwargs.keys(): setattr(self, key, kwargs[key]) 

  def _prerun(self):
    from os.path import join, exists
    from os import getcwd
    import shutil
    from . import files

#   for param in self:
#     print param.incar_string(self)
#   return

    # creates incar
    with open(join(self._tempdir, files.INCAR), "w") as incar:
      for param in self:
        incar.write(param.incar_string(self))
      incar.close()

    # creates kpoints file
    with open(join(self._tempdir, files.KPOINTS), "w") as kpoints:
      kpoints.write( self.kpoints(self) )
      kpoints.close()

    # creates poscar file
    print_poscar(self._system, self.species, self._tempdir)

    # creates POTCAR file
    with open(join(self._tempdir, files.POTCAR), 'w') as potcar:
      for s in self.__find_species__(_structure):
        cat, filename = "", ""
        if exists(join(s.path, "POTCAR")): cat, filename = "cat", join(s.path, "POTCAR")
        elif exists(join(s.path, "POTCAR.Z" )): cat, filename = "zcat", join(s.path, "POTCAR.Z")
        else: raise AssertionError, "Could not find potcar in " + s.path
        cmd = subprocess.Popen( [cat, os.path.join(s.path, filename) ], \
                                stdout = subprocess.PIPE )
        for line in cmd.stdout: print >>potcar, line[:len(line)-1]
      potcar.close()

    # checks for CHGCAR
    indir = self.indir
    if indir == None: indir = getcwd()
    with join(indir, files.CHGCAR) as chgcar:
      if exists(chgcar): shutil.copy(chgcar, self._tempdir)

    # checks for WAVECAR
    with join(indir, files.WAVECAR) as wavecar:
      if exists(chgcar): shutil.copy(chgcar, self._tempdir)

    # checks for EIGENVALUE
    with join(indir, files.EIGENVALUES) as eigenvalue:
      if exists(eigenvalue): shutil.copy(eigenvalue, self._tempdir)

  def _run(self):
     """ Isolates system call to vasp itself """
     from . import files
     with open(self.join(self._tempdir, files.STDOUT), "w") as stdout:
       with open(self.join(self._tempdir, files.STDERR), "w") as stderr:
         vasp_proc = subprocess.Popen( self.program, cwd = self._tempdir, \
                                       stdout = stdout, stderr = stderr, shell = True )
         vasp_proc.wait() # Wait for completion. Could set a timer here.

  def _postrun(self, repat, outdir):
     """ Copies files back to outdir """
     from os.path import exists, join, isdir
     from os import makedirs

     if not exists(outdir): makedirs(outdir)
     if not exists(outdir): raise IOError, "%s directory does not exist." % (outdir)
     if not isdir(outdir):  raise IOError, "%s is not a directory." % (outdir)

     notfound = []
     for filename in set(repat+files.minimal):
       if not exists(self._tempdir, filename):
         notfound.append(filename)
         continue
       copy( join(self._tempdir, filename), self.outdir )

     if len(notfound) != 0: raise IOError, "Files %s were not found.\n" % notfound

  def __call__(self, structure=None, outdir = None, repat = []):
    from ..opt.tempdir import Tempdir
    from ..opt.changedir import Changedir
    from os.path import exists, join
    from os import getcwd
    from shutil import copy2 as copy

    # set up
    if structure != None: self._system = structure
    elif not hasattr(self, "_system"): raise RuntimeError, "Internal bug.\n"
    if outdir == None: outdir = self.indir
    if outdir == None: outdir = getcwd()

    # creates temporary working directory
    workdir = self.workdir
    if workdir == None: workdir = getcwd()
    with Tempdir(workdir) as self._tempdir: 
      # creates INCAR and KPOINTS.
      # copies POSCAR, POTCAR, possibly WAVECAR and CHGCAR from indir
      self._prerun()

      # runs vasp.
      self._run()

      # now copies data back
      self._postrun(repat, outdir)

      # deletes system attribute.
      del self._system

