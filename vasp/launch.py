""" A class for single-shot VASP calculation """
from incar import Incar
from kpoints import Density


class Launch(Incar):
  """ A class to launch vasp """

  program = "vasp.mpi"
  """ shell command to launch vasp """

  POSCAR  = "POSCAR"
  """ Name of the input structure file """
  KPOINTS = "KPOINTS"
  """ Name of the kpoints file """
  INCAR   = "INCAR"
  """ Name of the input parameters file """
  POTCAR  = "POTCAR"
  """ Name of the pseudopotential file """
  WAVECAR = "WAVECAR"
  """ Name of the wavefunction file """
  CONTCAR = "CONTCAR"
  """ Name of the output structure file """
  CHGCAR  = "CHGCAR"
  """ Name of the output charge file """
  OSZICAR = "OSZICAR"
  """ Name of the energy minimization file """
  STDOUT  = "stdout"
  """ Name of the standard output file """
  STDERR  = "stderr"
  """ Name of the standard error file """
  EIGENVALUE = "EIGENVALUE"

  def __init__(self):
    Incar.__init__(self) 

    self.indir = ""
    self.workdir = ""
    self.kpoints = Density()

  def _prerun(self):
    from os.path import join, exists
    import shutil

#   for param in self:
#     print param.incar_string(self)
#   return

    # creates incar
    with open(join(self._tmpdir, self.INCAR), "w") as incar:
      for param in self:
        incar.write(param.incar_string(self))
      incar.close()

    # creates kpoints file
    with open(join(self._tmpdir, self.KPOINTS), "w") as kpoints:
      kpoints.write( self.kpoints(self) )
      kpoints.close()

    # creates poscar file
    print_poscar(self.system, self.species, self._tmpdir)

    # creates POTCAR file
    with open(join(self._tempdir, self.POTCAR) 'w') as potcar:
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
    with join(self.indir, self.CHGCAR) as chgcar:
      if exists(chgcar): shutil.copy(chgcar, self._tempdir)

    # checks for WAVECAR
    with join(self.indir, self.WAVECAR) as wavecar:
      if exists(chgcar): shutil.copy(chgcar, self._tempdir)

    # checks for EIGENVALUE
    with join(self.indir, self.EIGENVALUE) as eigenvalue:
      if exists(eigenvalue): shutil.copy(eigenvalue, self._tempdir)

  def __call__(self, structure, outdir = None, repat = []):
    from lada.opt import Tempdir, Changedir
    from os.path import exists, join
    from shutil import copy2 as copy

    # set up
    self.system = structure
    if outdir == None: outdir = self.indir

    # creates temporary working directory
    with Tempdir(self.workdir) as self._tempdir: 

      # creates INCAR and KPOINTS.
      # copies POSCAR, POTCAR, possibly WAVECAR and CHGCAR from indir
      self._prerun()

      # runs vasp.
      with open(self.join(self._tempdir, self.STDOUT), "w") as stdout:
        with open(self.join(self._tempdir, self.STDERR), "w") as stderr:
          vasp_proc = subprocess.Popen( self.program, cwd = self._tempdir, \
                                        stdout = stdout, stderr = stderr, shell = True )
          vasp_proc.wait() # Wait for completion. Could set a timer here.
          stderr.close()
        stdout.close()

      # now copies data back
      self._postrun(repat, outdir)

   def _postrun(self, repat, outdir)
      """ Copies files back to outdir """
      from os.path import exists, join, isdir
      from os import makedirs

      if not exists(outdir): makedirs(outdir)
      if not exists(outdir): raise IOError, "%s directory does not exist." % (outdir)
      if not isdir(outdir):  raise IOError, "%s is not a directory." % (outdir)

      repat.extend( [self.OUTCAR, self.INCAR, self.POSCAR, self.CONTCAR, \
                     self.STDOUT, self.STDERR, self.OSZICAR, self.EIGENVALUE ] )
      notfound = []
      for filename in set(repat):
        if not exists(self._tempdir, filename):
          notfound.append(filename)
          continue
        copy( join(self._tempdir, filename), self.outdir )

      if len(notfound) != 0: raise IOError, "Files %s were not found.\n" % notfound
