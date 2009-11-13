#
#  Version: $Id$
#

class Specie:
  """ Holds atomic specie information: 
      symbol is the atomic sumbol
      file is the path to the potcar directory.
           it will be assumed to be zipped.
  """

  def __init__(self, _symbol, _path ):
    """ Initializes a specie.
        _symbol is the atomic sumbol
        _path is the path to the potcar directory.
    """
    import os.path

    self.symbol = _symbol
    self.path = os.path.expanduser( _path )

class Vasp:
  """Can create 'vasp' jobs.
      Startup parameters:
      _ restart should be directory from which to restart, or nothing.
      _ species holds the list of Specie objects.
      _ iniwave can be either random or jellium. VASP manual recommends random.
      Electronic degrees of freedom:
      _ encut is the energy cutoff. if it is <= 0, then ENCUT = cutoff_safety * max(ENMAX)
        otherwise,  ENCUT = encut.
      _ cutoff_safety is a real number such that ECUT = cutoff_safety * EMAX,
        where EMAX is the maxium from the potcars of the species in the job.
      _ nelect is the number of electrons. 0 or "vasp" will have VASP compute it for us.
      _ nbands is the number of bands. 0 or "vasp" lest vasp figure it out.
      _ nspins is the number of spins, either 1 or 2.
      _ ng is the number of points NGX, NGY, NGZ.
        If it is zero, then runs vasp to find out best value.
        If it is negative, then do not use NG? tags.
      _ smearing is the type of smearing function used in metals. It can be
          "fermi x"  for Fermi-Dirac function.
          "gaussian x" for gaussians.
          "mp N x" for Methfessel-Paxton method, where N is the order of the expansion.
          "tretra x" for tetrahedron method without Blochl correction.
          "tretra bloech x" for tetrahedron method with Blochl correction.
          "metal x" is equivalent to "mp 1 0.2"
          "insulator" is equivalent to "tetra bloechl"
          x is the value of the smearing in eV. If it is not present, a value
          of 0.2eV is used by default.
       _ kpoint is the kpoint scheme. It can be:
          a list of atat.rVector3d or of 4-tuples with the direct coordinates and weights.
          "gamma N" to generate a gamma-centered NxNxN mesh
          "mp N " to generate an off-center Monkhorst-Pack NxNxN mesh.
          "file filename" to use a pre-generated file.
          "d=? o=?" to generate a mesh with density ? (from d=?) and an offset
                   (o=?). The algorithm is the same as vasp option 'auto',
                   except that an offset can be given. d is the number of
                   points per unit of reciprocal-space length.
          An empty string/array defaults to Gamma only.
      _ isym is the symmetrization method. It can be:
           off x,
           usp x recommended value for ultra-soft pseudos,
           paw x recommended value for paw pseudos
           default x will set isym depending on the pseudo of the species.
           x is an optional parameter specifying the precision.
      Electronic Minimization:
      _ algo can be "normal", "fast", or "very fast".
      _ ediff sets the convergence criteria of the charge SCF loop.
      _ prec sets the precision (cutoff, fft grids...). 
        It can be "low", "medium", "normal", "accurate", or "high".
      Ionic Minimization:
      _ relaxation specifies which degrees of freedom to relax: ionic, 
        cellshape, volume, and the method, which can be:
             md (IBRION = 0) molecular dynamics, no cell-shape optimization.
             local (IBRION = 1 ) rmm-diis is good close to the minimum
             global (IBRION = 2) conjugate-gradient is the most robust method.
      _ potim is the time-step for ionic-motions.
      _ nsw is the maximum number of ionic step. 40 by default.

      The path to the VASP executable should be specified in static variable Vasp.program.
  """

  program = "~/bin/vasp.mpi"
  files = ["POSCAR", "INCAR", "OUTCAR", "KPOINTS", "IBZKPT", \
           "CHGCAR", "WAVECAR", "CONTCAR", "DOSCAR", "OSZICAR", "POTCAR" ]

  def __init__(self):

    self.restart = ""
    self.species = [ Specie("Rb", "~/AtomicPotentials/pseudos/Rb_s/"), \
                     Specie("K", "~/AtomicPotentials/pseudos/K_s/") ]
    self.iniwave = "random"

    self.encut = -1
    self.cutoff_safety = 1.25
    self.nelect  = 0
    self.nbands = 0
    self.nspins = 1
    self.smearing = "gaussian 0.2"
    self.kpoint = "d=10 o=0.5"
    self.potim = 0.5
    self.isym  = "default 1e-5"

    self.algo = "fast"
    self.prec = "accurate"
    self.ediff = 1e-4
    self.fft = (-1,-1,-1)

    self.relaxation = "global volume ionic cellshape"
    self.nsw = 40

    # directory where infiles are stored.
    self.indir = ""

    self.other = {}

  def copy(self):
    """Returns a copy of this object
    """

    result = Vasp()
    result.restart = str(self.restart)
    result.species = list(self.species)
    result.iniwave = str(self.iniwave)
    result.encut = int(self.encut)
    result.fft = tuple(self.fft)
    result.cutoff_safety = float(self.cutoff_safety)
    result.nelect = float(self.nelect)
    result.nbands = int(self.nbands)
    result.nspins = int(self.nspins)
    result.smearing = str(self.smearing)
    result.kpoint = str(self.kpoint)
    result.potim = float(self.potim)
    result.isym = str(self.isym)
    result.algo = str(self.algo)
    result.prec = str(self.prec)
    result.ediff = float(self.ediff) 
    result.relaxation = str(self.relaxation)
    result.nsw = int(self.nsw)
    result.indir = str(self.indir)
    return result


  def read(self, _filename):
    """Attempts ro recreate a Vasp oject from path _filename, using INCAR or OUTCAR.
    """
    import os

    if os.path.isfile( _filename ):
      if os.path.basename( _filename) == "INCAR": return self.read_incar( _filename )
      elif os.path.basename( _filename) == "OUTCAR": return self.read_outcar( _filename )
      else: raise RuntimeError, "Unknown file type " + _filename + "\n" 
    if os.path.exists( os.path.join(_filename, "OUTCAR") ):
      return self.read_outcar( os.path.join(_filename, "OUTCAR") )
    elif os.path.exists( os.path.join(_filename, "INCAR") ):
      return self.read_incar( os.path.join(_filename) )
    else: raise RuntimeError,   "Could not find INCAR or OUTCAR in directory " \
                              + _filename + "\n" 

  def read_incar( self, _filename ):
    """Attempts ro recreate a Vasp oject from path _filename, using INCAR.
    """

    result = Vasp()
    filename = str(_filename)
    if os.path.isdir( _filename ):
      filename  = os.path.join(_filename, "INCAR")
      if not os.path.exists( filename ):
        raise RuntimeError, "Could not find INCAR in " + _filename + ".\n"
    result.indir = str( os.path.split(filename)[0] )
    result.species = list(self.species)
    file.open( filename, 'w' )
    encut = re.compile( "ENCUT\s*=\s*(\S+)" )
    ispin = re.compile( "ISPIN\s*=\s*(\S+)" )
    ispin = re.compile( "ISPIN\s*=\s*(\S+)" )
    for line in file:
      p = encut.search( line ) 
      if p != None:
        result.encut = float( p.group(1) )
        continue
      p = ispin.search( line ) 
      if p != None:
        result.ispin = int( p.group(1) )
        continue
      p = system.search( line ) 
      if p != None:
        result.ispin = int( p.group(1) )
        continue



  def potcar(self, _structure, _filename ):
    """Returns a string containing the POTCAR.
    """
    import os.path
    import subprocess

    file = open( _filename, 'w')
    for s in self.__find_species__(_structure):
      cat = ""
      if os.path.exists( os.path.join(s.path, "POTCAR") ):  
        p = open( os.path.join(s.path, "POTCAR"), 'r' )
        for line in os.path.join(s.path, "POTCAR"): print >>file, line[:len(line)-1]
        p.close()
      elif os.path.exists( os.path.join(s.path, "POTCAR.Z" )): 
        cmd = subprocess.Popen( ["zcat", os.path.join(s.path, "POTCAR") ], \
                                stdout = subprocess.PIPE )
        for line in cmd.stdout: print >>file, line[:len(line)-1]
      else: raise AssertionError, "Could not find potcar in " + s.path
    file.close()

  def kpoints(self, _structure):
    """Returns a string with the KPOINTS file.
    """
    import re
    import os.path
    from lada import atat
    from math import sqrt
    
    if len(self.kpoint) == 0: 
      return "\n1\nGamma\n1 1 1\n0 0 0\n"
    elif isinstance(self.kpoint, str):
      if re.search("g(amma)?", self.kpoint.lower()):
        N = 1
        if len( self.kpoint.split() ) > 1:
          N = int( self.kpoint.split()[1] )
        return "\n0\nGamma\n%i %i %i\n0 0 0\n" % ( N, N, N )
      elif re.search("m(p)?", self.kpoint.lower()):
        N = 1
        if len( self.kpoint.split() ) > 1:
          N = int( self.kpoint.split()[1] )
        return "\n0\nM\n%i %i %i\n0 0 0\n" % ( N, N, N )
      elif re.search("d\s*=\s*?(\d+)", self.kpoint): 
        N = re.search("d\s*=\s*(\d+)", self.kpoint)
        N = int( N.group(1) )
        offset = re.search("o\s*=\s*(\S+)", self.kpoint) 
        if offset == None: 
          return "\n0\nAuto\n%i\n" % (N)
        else:
          offset = float(offset.group(1))
          recip = atat.transpose( atat.inverse(_structure.cell) ) / _structure.scale
          string = "\n0\nM\n"
          for i in range(0,3):
            b = sqrt( atat.norm2( atat.rVector3d( recip[(0,i)], recip[(1,i)], recip[(2,i)]) ) )
            string += "%i " % ( int(max( 1, N*b+0.5)) )
          string += "\n%f %f %f" % (offset, offset, offset)
          return string
      elif re.search("max\s*=\s*?(\d+)", self.kpoint): 
        N = re.search("d\s*=\s*(\d+)", self.kpoint)
        N = int( N.group(1) )
        recip = atat.transpose( atat.inverse(_structure.cell) ) / _structure.scale
        string = "\n0\nM\n"
        m = 0
        for i in range(0,3):
          b = sqrt( atat.norm2( atat.rVector3d( recip[(0,i)], recip[(1,i)], recip[(2,i)]) ) )
          b = int(max( 1, N*b+0.5))
          if m < b: m = b;
        string += "%i %i %i\n0 0 0" % (m, m, m)
        return string
      elif re.search("file (\S+)", self.kpoint.lower()):
        filename = self.kpoint().split()[1] 
        if not os.path.exists(filename):
          raise RuntimeError, "Could not find file " + filename + ".\n"
        file = open(filename, 'r')
        result = ""
        for line in file:
          result += line
        file.close()
        return result
    else:
      result = "\n%i\n" % len(self.kpoint)
      for k in self.kpoint:
        if len(k) == 2:
          result += "%s %f\n" % (atat.rVector3d(k[0]), k[1])
        elif len(k) == 4:
          result += "%s %f\n" % (atat.rVector3d(k[:3]), k[1])
        else: raise RuntimeError, "Uknown kpoint format %s\n" % (k)
      return result


  def figure_fft( self, _structure, _header, _footer ):

    import tempfile 
    import os
    import shutil
    import subprocess
    import re

    vasp = self.copy()
    vasp.indir = tempfile.mkdtemp()
    vasp.fft = (-1, -1, -1)
    vasp.other["NELMDL"] = 0
    vasp.other["NELM"] = 0
    vasp.other["NELMIN"] = 0
    vasp.relaxation = "static"
    script = vasp.prepare(_structure, wpath=vasp.indir, header=_header, footer=_footer, repat=[] )
    filename = os.path.join(vasp.indir, 'script')
    file = open( filename, 'w' )
    print >>file, script
    file.close()
    filename = os.path.join(vasp.indir, 'script')
    call = subprocess.Popen(["/bin/sh", filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in call.stdout:
      print line[:(len(line)-1)]
    if not os.path.exists( os.path.join(vasp.indir, "OUTCAR") ):
      print vasp.indir
      raise RuntimeError, "Could not figure NG?.\n" 
    file = open( os.path.join(vasp.indir, "OUTCAR"), 'r' )
    for line in file:
      if re.search("I would recommend the setting", line): break;
    fft = [0,0,0]
    ng_regex = re.compile("WARNING: wrap around error must be expected set NG(X|Y|Z) to\s+(\d+)")
    g_regex = re.compile("NG(X|Y|Z) is ok and might be reduce to\s+(\d+)")
    allset = 0
    multiple = 8
    for line in file:
      p = ng_regex.search(line)
      if p != None:
        if p.group(1) == 'X':
	  fft[0] = int(p.group(2)) 
          if fft[0] % multiple: fft[0] += multiple - fft[0] % multiple
          allset += 1
        elif p.group(1) == 'Y':
          fft[1] = int(p.group(2))
          if fft[1] % multiple: fft[1] += multiple - fft[1] % multiple
          allset += 1
        elif p.group(1) == 'Z':
          fft[2] = int(p.group(2))
          if fft[2] % multiple: fft[2] += multiple - fft[2] % multiple
          allset += 1
        if allset == 3: break;
        continue;
      p = g_regex.search(line)
      if p != None:
        if p.group(1) == 'X':
          fft[0] = int(p.group(2)) 
          if fft[0] % multiple: fft[0] += multiple - fft[0] % multiple
          allset += 1
        elif p.group(1) == 'Y':
          fft[1] = int(p.group(2))
          if fft[1] % multiple: fft[1] += multiple - fft[1] % multiple
          allset += 1
        elif p.group(1) == 'Z':
          fft[2] = int(p.group(2))
          if fft[2] % multiple: fft[2] += multiple - fft[2] % multiple
          allset += 1
        if allset == 3: break;


    file.close()
    shutil.rmtree( vasp.indir )
    return fft
    



  def incar(self, _structure):
    """Returns a string with the INCAR
    """
    import os.path
    import re
    from lada import atat, physics
    from math import sqrt
    # checks existence of pseudos.
    species = self.__find_species__(_structure)

    # how and where to start
    istart = "0   # start from scratch"
    icharg = "2   # superpositions of atomic densities"
    if len(self.restart) == 0:
      istart = "0   # start from scratch"
    elif not os.path.exists( self.restart ):
      raise RuntimeError, "Could not find restart directory " + self.restart + "\n";
    else:
      ewave = os.path.exists( os.path.join( self.restart, 'WAVECAR' ) )
      echarge = os.path.exists( os.path.join( self.restart, 'CHGCAR' ) )
      if ewave:
        istart = "1  # restart"
        icharg = "0   # from wavefunctions " + os.path.join( self.restart, 'WAVECAR' )
      elif echarge:
        istart = "1  # restart"
        icharg = "1   # from charge " + os.path.join( self.restart, 'CHGCAR' )
      else: 
        istart = "0   # start from scratch"
        icharg = "2   # superpositions of atomic densities"
    iniwave = 0
    if  istart[0] == '0':
      if self.iniwave.lower() == "random": 
        iniwave = "1   # starts from random waves."
      elif self.iniwave.lower() == "jellium":
        iniwave = "0   # starts from jellium waves."
      else: raise RuntimeError, "Uknown iniwave type: " + self.iniwave + "\n" 

    result = "SYSTEM = %s\n\n" % (_structure.name)
    result += "# how to start\n"
    result += "ISTART = %s\n" % (istart)
    result += "ICHARG = %s\n" % (icharg)
    if istart[0] == '0': result += "INIWAV = %s\n\n" % (iniwave)


    encut = 0
    if self.encut <= 0: 
      a = self.__get_encut__( self.species )
      encut =    " %5.2f     # cutoff: %f * %f" \
               % ( self.cutoff_safety * a, self.cutoff_safety, a )
    else: encut = "%f" % (self.encut)
    algo = "fast"
    if self.algo.lower() == "normal": algo = "Normal     # Minization algorithm"
    elif self.algo.lower() == "fast": algo = "Fast       # Minization algorithm"
    elif self.algo.lower().replace("_", " ") == "very fast":
      algo = "Very_Fast # Minization algorithm"
    else: raise RuntimeError, "Uknown algorithm: " + self.algo + "\n" 
    if self.nspins != 1 and self.nspins != 2: 
      raise RuntimeError, "No exotic physics please. spin = %i?" % (self.nspins)



    prec = "accurate"
    if self.prec.lower() == "low":
      prec = "Low        # energy cutoff, fft grid points, real-space points.."
    elif self.prec.lower() == "medium":
      prec = "Medium     # energy cutoff, fft grid points, real-space points.."
    elif self.prec.lower() == "accurate":
      prec = "Accurate   # energy cutoff, fft grid points, real-space points.."
    elif self.prec.lower() == "high":
      prec = "High       # energy cutoff, fft grid points, real-space points.."
    else: raise RuntimeError, "Uknown precision: " + self.prec + "\n" 

    smearing = ""
    sigma = 0.2
    if re.search("fermi", self.smearing.lower()):
      smearing = "-1         # Fermi-Dirac smearing."
      if len( self.smearing.split() ) > 1:
        sigma = float(self.smearing.split()[1])
    elif re.search("gaussian", self.smearing.lower()):
      smearing = "0          # Gaussian smearing."
      if len( self.smearing.split() ) > 1:
        sigma = float(self.smearing.split()[1])
    elif re.search("mp", self.smearing.lower()):
      N = 1
      if len( self.smearing.split() ) > 1:
        N = int(self.smearing.split()[1])
        if len( self.smearing.split() ) > 2:
          delta = float(self.smearing.split()[2])
      smearing =   "%2i          # Methfessel-Paxtion smearing of order %i."\
                 % (N, N)
    elif re.search("bloechl", self.smearing.lower()):
      smearing = "-5         # Tetrahedron method with Bloechl correction."
      if len( self.smearing.split() ) > 1:
        sigma = float(self.smearing.split()[1])
    elif re.search("tetra", self.smearing.lower()):
      smearing = "-4         # Tetrahedron method without Bloechl correction."
      if len( self.smearing.split() ) > 1:
        sigma = float(self.smearing.split()[1])
    elif re.search("insulator", self.smearing.lower()):
      smearing = "-5         # Tetrahedron method with Bloechl correction."
      if len( self.smearing.split() ) > 1:
        sigma = float(self.smearing.split()[1])
    elif re.search("metal", self.smearing.lower()):
      smearing = "1          # Tetrahedron method with Bloechl correction."
      if len( self.smearing.split() ) > 1:
        sigma = float(self.smearing.split()[1])
    else: raise RuntimeError, "Unkown smearing " + self.smearing + "\.n"

    isym = ""
    if re.search( "off", self.isym.lower() ):
      isym = "0          # No symmetrization"
    elif re.search( "ups", self.isym.lower() ):
      isym = "1          # Symmetrization for UPS."
    elif re.search( "paw", self.isym.lower() ):
      isym = "3          # Symmetrization for UPS."
    elif re.search( "default", self.isym.lower() ):
      isym = "default"
    else: raise RuntimeError, "Unkown symmetrization " + self.isym + "\.n"
    symprec = 1e-5
    if len(self.isym.split()) > 1:
      symprec = float( self.isym.split()[1])



    result += "# Electronic degrees of freedom.\n"
    result += "ENCUT   = %s\n" % (encut)
    if self.fft[0] > 0: result += "NGX     = %s\nNGY     = %s\nNGZ     = %s\n" % tuple(self.fft)
    result += "ISPIN   = %i\n" % (self.nspins)
    if int(self.nelect) > 0:
      result += "NELECT  = %3i        "\
                "# Number of electrons (e.g. charge)\n" % (int(self.nelect))
    if int(self.nbands) > 0:
      result += "NBANDS  = %3i        "\
                "  # Number of bands (e.g. not charge)\n" % (int(self.nbands))
    result += "ISMEAR  = %s\nSIGMA   = %f\n" % (smearing, sigma)
    if isym != "default":
      result += "ISYM    = %s\n" % (isym)
    result += "SYMPREC = %3.2e   # Precision to determine symmetry operation.\n" % (symprec)
    if int(self.nelect) <= 0: result += "  # using default number of electrons from VASP.\n"
    if int(self.nbands) <= 0: result += "  # using default number of bands from VASP.\n"
    if isym == "default":     result += "  # using default ISYM from VASP\n"

    result += "\n# Electronic minimization.\n"
    result += "ALGO    = %s\n" % (algo)
    result +=   "EDIFF   = %3.2e   # convergence criteria for the SCF loop\n"\
              % (self.ediff / float(len(_structure.atoms)) )
    result += "PREC    = %s\n" % (prec)



    isif = 0
    ionic =  re.search( "ionic", self.relaxation.lower() ) != None
    cellshape = re.search( "cell(\s+|-|_)?shape", self.relaxation.lower() ) != None
    volume = re.search( "volume", self.relaxation.lower() ) != None
    if (not ionic) and (not cellshape) and (not volume) :
      isif = "0                # static calculation " 
    elif ionic and (not cellshape) and (not volume) :
      isif = "2                # relaxing atomic positions only."
    elif ionic and cellshape and volume: 
      isif = "3                # relaxing all structural degrees of freedom."
    elif ionic and cellshape and (not volume):
      isif = "4                # relaxing atomic positions and cell-shape at constant volume."
    elif(not ionic) and  cellshape and (not volume):
      isif = "5                # relaxing cell-shape at constant atomic-positions and volume."
    elif(not ionic) and  cellshape and volume:
      isif = "6                # relaxing volume and cell-shape at constant atomic-positions."
    elif(not ionic) and (not cellshape) and volume:
      isif = "7                # relaxing volume only."
    elif ionic and (not cellshape) and volume:
      raise RuntimeError, "VASP does not allow relaxation of atomic position"\
                          "and volume at constant cell-shape.\n"
    ibrion = ""
    if len( self.relaxation ) == 0 or self.relaxation == "static":
      ibrion = "0                # static calculation."
    elif re.search("md", self.relaxation.lower()):
      if cellshape or volume:
        raise RuntimeError, "Cannot perform cell-shape/volume relaxation " \
                            "during molecular dynamics.\n"
      ibrion = "0                # Molecular dynamics"
    elif re.search("local", self.relaxation.lower()):
      ibrion = "1                # Quasi-Newton, best when close to optimum."
    elif re.search("global", self.relaxation.lower()):
      ibrion = "2                # Conjugate-gradient, robust when far from optimum."

    result += "\n# Ionic minimization.\n"
    result += "IBRION = %s \n"  % ( ibrion )
    result += "ISIF   = %s \n"  % ( isif )
    result += "POTIM  = %4.3e        "\
              "# \"time\"-scale of ionic relaxation.\n" % (self.potim)
    if ibrion[0] != '0' and ibrion[:2] != "-1":
      result += "NSW    = %4i          "\
                "# Number of ionic steps.\n" % (self.nsw)

    if len(self.other) != 0:
      result += "\n# Other tags.\n"
      for key in self.other.keys():
        result += "%s = %s\n" % ( key, self.other[key] )

    return result

  def __find_species__( self, _structure ):
    """Returns a list of species in the structure.
    """
    import os.path

    results = []
    for atom in _structure.atoms:
      for s in self.species:
        if s.symbol == atom.type and not(s in results):
          results.append( s )
          a = os.path.join( s.path, "POTCAR" )
          b = os.path.join( s.path, "POTCAR.Z" )
          if not (os.path.exists(a) or os.path.exists(b)):
            raise AssertionError, "Could not find potcar in " + s.path
    return results


  def __get_encut__( self, _species ):
    """ Retrieves max ENMAX from list of species.
    """
    import os.path
    import subprocess
    import re
    from math import ceil

    result = 0
    for s in _species: 
      stdout = ""
      if os.path.exists( os.path.join(s.path, "POTCAR") ):
        filename = os.path.join(s.path, "POTCAR")
        cmd = subprocess.Popen(["grep", "ENMAX", s.path + "POTCAR"], \
                               stdout=subprocess.PIPE)
        stdout = cmd.stdout.readline()
      elif os.path.exists( os.path.join(s.path, "POTCAR.Z") ):
        filename = os.path.join(s.path, "POTCAR.Z")
        cmd0 = subprocess.Popen(["zcat", filename], \
                                stdout=subprocess.PIPE)
        cmd = subprocess.Popen(["grep", "ENMAX"], \
                               stdin=cmd0.stdout, stdout=subprocess.PIPE)
        stdout = cmd.communicate()[0]
      else: raise AssertionError, "Could not find potcar in " + s.path
  
      r = re.compile("ENMAX\s+=\s+(\S+);\s+ENMIN")
      p = r.search(stdout)
      if p == None: raise AssertionError, "Could not retrieve ENMAX from " + s.path
      if result < float( p.group(1) ): result = float( p.group(1) )

    return ceil(result)
    
  def prepare(self, _structure, wpath = "", repat = "all", header = "", footer="" ):
    """Prepares a script for launching vasp.
       _structure is the structure over which to perform calculations.
       _repat are the file to repatriate: "outcar, poscar, chgcar, wavecar, contcar, ibzkpt"
              "all" will copy everything back to _outpath.
              "none" will do the opposite
       _wpath is the directory where to compute the structure.
    """
    import os
    import tempfile
    from lada import crystal


    if len( self.indir ) == 0:
      if len(_structure.name) != 0:
        tempfile.tempdir = os.getcwd()
        self.indir = tempfile.mkdtemp(prefix=_structure.name + "_")
      else:
        tempfile.tempdir = os.getcwd()
        self.indir = tempfile.mkdtemp()
    elif not os.path.exists(self.indir): os.mkdir( self.indir )

    if len(wpath) == 0: wpath = str(self.indir)
    if not os.path.exists(wpath): os.mkdir(wpath)
    tempfile.tempdir = wpath

    delfft = 0
    if self.fft[0] == 0:
      self.fft = self.figure_fft( _structure, header, footer ) 
      delfft = 1


    # create INCAR
    file = open( os.path.join(self.indir, "INCAR"), 'w')
    print >>file, self.incar( _structure )
    file.close()

    if delfft: self.fft = (0,0,0)
    # create KPOINTS
    file = open( os.path.join(self.indir, "KPOINTS"), 'w')
    print >>file, self.kpoints( _structure )
    file.close()
    # create POTCAR
    self.potcar( _structure, os.path.join(self.indir, "POTCAR") )
    # create POSCAR
    crystal.print_poscar\
    (\
      _structure, 
      tuple( s.symbol for s in self.__find_species__(_structure) ), 
      self.indir
    )
    
    result  = str(header) + "\n"
    # copy files from indir.
    if os.path.abspath( self.indir ) != os.path.abspath( wpath ):
      result += "\n# copy vasp files\n"
      for file in Vasp.files:
        if os.path.exists( os.path.join(self.indir, file) ):
          result += "cp " + os.path.join( self.indir, file ) + " " + wpath + "\n";
    # copy charge and/or wave file from restart directory.
    if len(self.restart) != 0:
      if     ( os.path.abspath( self.restart ) != os.path.abspath( wpath ) )\
         and ( os.path.abspath( self.restart ) != os.path.abspath( self.indir ) ):
        result += "\n# copy restart files\n"
        for file in ['CHGCAR', 'WAVECAR']:
          if os.path.exists( os.path.join(self.restart, file) ):
            result += "cp " + os.path.join( self.restart, file ) + " " + wpath + "\n";
  
    # run program
    result += "\n# Go to execution path and start run\ncd " + wpath + "\n"
    result += Vasp.program + "\n"

    # copy files back to indir.
    if os.path.abspath( self.indir ) != os.path.abspath( wpath ):
      if isinstance(repat, str):
        if repat.lower() == "all":
          files = Vasp.files
        else: files = repat.split()
      else: files = repat
      for file in files:
        result += "cp " + os.path.join(wpath, file) + " " + os.path.abspath(self.indir) + "\n\n";

    result += footer + "\n";
    
    return result


  def final(self):
    """ Returns structure at end of run.
    """
    import os.path
    from lada import crystal
    import re

    filename = os.path.join( self.indir, "CONTCAR" )
    if not os.path.exists( filename ):
      raise RuntimeError, "Could not find " + filename
    filename = os.path.join( self.indir, "OSZICAR" )
    if not os.path.exists( filename ): 
      raise RuntimeError, "Could not find " + filename
    filename = os.path.join( self.indir, "OUTCAR" )
    if not os.path.exists( filename ): 
      raise RuntimeError, "Could not find " + filename

    # find the species.
    file = open( filename, 'r')
    aset = set()
    for line in file:
      if re.search("POTCAR:", line) != None: break
    aset.add( line.split()[2] )
    for line in file:
      if re.search("POTCAR:", line) == None: break
      aset.add( line.split()[2] )
    file.close()
 
    species = []
    for s in self.species:
      if s.symbol in aset: species.append( s.symbol ) 
    for s in species:
      if not ( s in [ u.symbol for u in self.species ] ):
        raise RuntimeError, "Unknown specie " + s + " found in OUTCAR.\n"
 
    # reads (last) structure
    result = crystal.sStructure()
    filename = os.path.join( self.indir, "CONTCAR" )
    result = crystal.read_poscar( tuple(species), filename )
 
    # reads (last) energy.
    filename = os.path.join( self.indir, "OSZICAR" )
    file = open( filename, 'r')
    found = 0
    for line in file:
      if re.search("E0", line) == None: continue;
      found = 1
      result.energy = float(line.split()[4])
    if not found:
      raise RuntimeError, "Could not find energy.\n";

    return result



def main():
  from lada import atat, crystal
  import os.path
  import shutil
  import subprocess


  structure = crystal.sStructure()
  structure.cell = atat.rMatrix3d( [[1,0,0],[0,1,0],[0,0,1]] )
  atom = crystal.StrAtom()
  atom.type = "Rb"
  atom.pos = atat.rVector3d([0,0,0])
  structure.atoms.append( crystal.StrAtom(atom) )
  atom.type = "K"
  atom.pos = atat.rVector3d([0.5,0.5,0.5])
  structure.atoms.append( crystal.StrAtom(atom) )
  structure.name = "KRb"
  structure.scale = 6

  K = Specie( "K", "~/AtomicPotentials/pseudos/K_s" )
  Rb = Specie( "Rb", "~/AtomicPotentials/pseudos/Rb_s" )

  vasp = Vasp()
  vasp.indir = "KRb"
  if os.path.exists( vasp.indir ):  shutil.rmtree( vasp.indir )
  script = vasp.prepare( structure, wpath="KRb", header="module load mpigm.pgi" )
  file = open( os.path.join( vasp.indir, 'script' ), 'w' )
  print >>file, script
  file.close()
  subprocess.call( ["bash", os.path.join(vasp.indir, "script") ] )

# structure = vasp.final()
# print structure.energy, structure.name
# print structure


if __name__ == "__main__":
  main()
