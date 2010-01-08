from incar import Incar
from kpoints import Density


class Run(Incar):
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

      repat.extend( [self.OUTCAR, self.INCAR, self.POSCAR, self.CONTCAR, \
                     self.STDOUT, self.STDERR, self.OSZICAR, self.EIGENVALUE ] )
      notfound = []
      for filename in set(repat):
        if not exists(self._tempdir, filename):
          notfound.append(filename)
          continue
        copy( join(self._tempdir, filename), self.outdir )

      assert len(notfound) == 0, "Files %s were not found.\n" % notfound

class Kpoint(object):
  """ A structure defining eigenvalues a point in reciprocal space 
  
      A kpoint  consists of cartesian coordinates, a weight (multiplicity), the
      eigenvalues and occupations. 
  """
  
  def __init__(self, pos, weight, nbands)
    from numpy import vector
    self.pos = pos
    self.weight = weight
    self.eigenvalues = vector([0e0 for u in range(nbands)], dtype="float64")

class Extract(object)
  """ Extracts VASP values """

  def __init__(self, indir = ""): self.indir = indir

  def _get_energy_sigma0(self):
    """ Gets total energy extrapolated to $\sigma=0$ from vasp run """
    from os.path import exists, join
    import re

    path = Run.OUTCAR 
    if self.indir != "": path = join(self.indir, path)
    assert exists(path), "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = re.compile( r"""energy\s+without\s+entropy\s*=\s*(\S+)\s+
                              energy(sigma->0)\s+=\s+(\S+)""", re.X)
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(2))
    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  energy_sigma0 = property(_get_energy_sigma0)
  r""" Gets total energy extrapolated to $\sigma=0$ from vasp run """

  def _get_energy(self):
    """ Gets total energy from vasp run """
    from os.path import exists, join
    import re

    path = Run.OUTCAR 
    if self.indir != "": path = join(self.indir, path)
    assert exists(path), "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = re.compile( r"""energy\s+without\s+entropy\s*=\s*(\S+)\s+
                              energy(sigma->0)\s+=\s+(\S+)""", re.X)
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))
    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  energy = property(_get_energy)
  r""" Gets total energy from vasp run """

  def _get_free_energy(self):
    """ Gets total free energy from vasp run """
    from os.path import exists, join
    import re

    path = Run.OUTCAR 
    if self.indir != "": path = join(self.indir, path)
    assert exists(path), "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = re.compile( r"""free\s+energy\s+TOTEN\s*=\s*(\S+)\s+eV""" )
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))

    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  free_energy = property(_get_free_energy)
  r""" Gets total free energy from vasp run """

  def _get_fermi_energy(self):
    """ Gets total free energy from vasp run """
    from os.path import exists, join
    import re

    path = Run.OUTCAR 
    if self.indir != "": path = join(self.indir, path)
    assert exists(path), "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = re.compile( r"""E-fermi\s*:\s*(\S+)""" )
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))

    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  fermi_energy = property(_get_fermi_energy)
  r""" Gets fermi energy from vasp run """

  def _get_structure(self):
    """ Gets structure from CONTCAR file and total energy from OUTCAR """
    from os.path import exists, join
    import re

    path = Run.CONTCAR 
    if self.indir != "": path = join(self.indir, path)
    assert exists(path), "File %s does not exist.\n" % (path)

    result = crystal.read_poscar(tuple(self.species), self.CONTCAR)
    result.energy = self.energy

    return result
  system = property(_get_structure)
  r""" Fills a crystal.sStructure from CONTCAR and free energy in OUTCAR. """

  def _get_fft(self):
    """ Returns recommended or actual fft setting """
    from os.path import exists, join
    import re

    path = Run.OUTCAR 
    if self.indir != "": path = join(self.indir, path)
    assert exists(path), "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:

      # find start
      for line in file:
        if re.search("I would recommend the setting", line): break;
      ng_regex = re.compile(r"WARNING:\s+wrap\s+around\s+error\s+
                              must\s+be\s+expected\s+set\s+NG(X|Y|Z)\s+to\s+(\d+)", re.X)
      g_regex = re.compile(r"NG(X|Y|Z)\s+is\s+ok\s+and\s+might\s+be\s+reduce\s+to\s+(\d+)", re.X)
      found_regex = re.compile(r"""dimension\s+x,y,z\s+
                                   NGX\s+=\s+(\d+)\s+
                                   NGY\s+=\s+(\d+)\s+
                                   NGZ\s+=\s+(\d+)""", re.X)

      allset = 0
      fft = None, None, None
      for line in file:
        p = ng_regex.search(line)
        if p != None:
          if p.group(1) == 'X':
            fft[0] = int(p.group(2)) 
            allset += 1
          elif p.group(1) == 'Y':
            fft[1] = int(p.group(2))
            allset += 1
          elif p.group(1) == 'Z':
            fft[2] = int(p.group(2))
            allset += 1
          if allset == 3: break;
          continue;
        p = g_regex.search(line)
        if p != None:
          if p.group(1) == 'X':
            fft[0] = int(p.group(2)) 
            allset += 1
          elif p.group(1) == 'Y':
            fft[1] = int(p.group(2))
            allset += 1
          elif p.group(1) == 'Z':
            fft[2] = int(p.group(2))
            allset += 1
          if allset == 3: break;
          continue
        p = found_regex.search(line)
        if p != None:
          fft = ( int(p.group(1)), int(p.group(2)), int(p.group(3)) )
          break;

      assert fft[0] != None, "File %s is incomplete or incoherent.\n" % (path)
      assert fft[1] != None, "File %s is incomplete or incoherent.\n" % (path)
      assert fft[2] != None, "File %s is incomplete or incoherent.\n" % (path)

      multiple = 8
      for i in range(3):
        if fft[i] % multiple: fft[i] += multiple - fft[i] % multiple
      return fft
    raise RuntimeError, "File %s could not be opened.\n" % (path)
  fft = property(_get_fft)
  r""" Gets recommended fft grid (for wavefunctions). """


      

# def _get_eigenvalues(self):
#   """ Gets eigenvalues from OUTCAR.

#       Returns a list of Kpoints objects. 
#   """
#   from os.path import exists, join
#   import re

#   path = Run.OUTCAR 
#   if self.indir != "": path = join(self.indir, path)
#   assert exists(path), "File %s does not exist.\n" % (path)
#   
#   with open(path, "r") as file:
#     # first gets number of kpoints and bands.
#     regex = re.compile( r"""k-Points\s+NKPTS\s*=\s*(\d+)\s*
#                             number\s+of\s+bands\s+NBANDS\s+=\s+(\d+)""", re.X )
#     nkpt, nbands = None, None
#     for line in file:
#       match = regex.search(line)
#       if match != None: 
#         nkpt, nbands = float(match.group(1)), float(match.group(2))
#         break
#     if nkpt == None: raise RuntimeError, "File %s is incomplete.\n" % (path)

#     # now reads kpoints 
#     results = []
#     file.seek(0)
#     regex = re.compile( r"""k-point\s+(\d+)\s+:\s+(\S+)\s+(\S+)\s+(\S+)""")
#     eigregex = re.compile( r"""\s+(\d+)\s+(\S+)\s+(\S+)""" )
#     last_read = -1 
#     is_reading = None
#     kpoint = None
#     for line in file:
#       if is_reading != None
#         if is_reading == 0: is_reading +=1 
#         elif is_reading <= nbands:
#           match = eigregex.search(line)
#           assert match != None, "File %s is incomplete or incoherent" % (path)
#           assert int(match.group(1)) == is_reading, "File %s is incomplete or incoherent" % (path)
#           result[kpoint].eigenvalues[is_reading-1, 0] =  float(match.group(2))
#           result[kpoint].eigenvalues[is_reading-1, 1] =  float(match.group(3))
#         else: 
#           assert is_reading == nbands + 1, "File %s is incomplete or incoherent" % (path) 
#           is_reading = None
#           last_read = kpoint
#           continue
#         is_reading +=1
#       else:
#         match = regex.search(line)
#         if match == None: continue
#         kpoint = int(match.group(1)) - 1
#         is_reading = 0
#         if last_read == nkpts-1: last_read = -1
#         assert last_read + 1 == kpoint, "File %s is incomplete or incoherent" % (path) 
            






      



if __name__ == "__main__":
  from lada import crystal, atat
  from specie import Specie
  
  vasp = Vasp() 
  vasp.species = [Specie("Al", "~/AtomicPotentials/pseudos/K_s")]
  vasp.fftgrid.value = (10,10,10)
  structure = crystal.sStructure()
  structure.scale = 1e0
  structure.cell = atat.rMatrix3d([[2,0.5,0.5],[0.0,0,0.5],[0.0,0.5,0]])
  structure.atoms.append( crystal.StrAtom() )
  structure.atoms[0].pos = atat.rVector3d(0,0,0)
  structure.atoms[0].type = "Al"

  print structure

  vasp(structure)







