""" Subpackage containing extraction methods for CRYSTAL output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract', 'MassExtract']

from ..tools import make_cached
from ..tools.extract import search_factory, AbstractExtractBase
from ..jobfolder import AbstractMassExtract
from ..error import GrepError
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'crystal.out')

class ExtractBase(object):
  """ Implementation class for extracting data from CRYSTAL output """

  def __init__(self):
    """ Initializes the extraction class. """
    super(ExtractBase, self).__init__()

  @property 
  @make_cached
  def crystal_family(self):
    """ Crystal family """
    result = self._find_first_STDOUT(r"^\s*CRYSTAL FAMILY\s*:\s*(\S+)")
    if result is None:
      raise GrepError( 'Could not grep crystal family from '                   \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    return result.group(1).lower()
  @property 
  @make_cached
  def crystal_class(self):
    """ Crystal class """
    result = self._find_first_STDOUT(r"^\s*CRYSTAL CLASS\s*(?:.+):\s*(.+)$")
    if result is None: 
      raise GrepError( 'Could not grep crystal class from'                     \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    return result.group(1).lower().rstrip().lstrip()

  @property 
  @make_cached
  def centrosymmetric(self):
    """ Crystal class """
    result = self._find_first_STDOUT(r"^\s*SPACE GROUP\s*\((\S+)\)\s*:")
    if result is None:
      raise GrepError( 'Could not grep centro-symmetricity from '              \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    return result.group(1).lower().rstrip().lstrip() != 'noncentrosymmetric'

  @property 
  @make_cached
  def space_group(self):
    """ Crystal class """
    result = self._find_first_STDOUT(r"^\s*SPACE GROUP\s*(?:.+)\s*:\s*(.+)")
    if result is None:
      raise GrepError( 'Could not grep space-group from '                      \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    return result.group(1).rstrip().lstrip()

  @property
  @make_cached
  def direct_lattice_cell(self):
    """ Direct lattice cell. """
    from numpy import array
    try: 
      file = self.__stdout__()
      for line in file:
        if len(line) <= 55: continue
        if line.split() == [ 'DIRECT', 'LATTICE', 'VECTORS', 'CARTESIAN',
                             'COMPONENTS', '(ANGSTROM)']: break
      file.next()
      return array( [file.next().split() for u in xrange(3)],                  \
                    dtype='float64' ).T
    except StopIteration:
      raise GrepError( 'Reached end of file whend grepping '                   \
                       'direct lattice cell vectors in '                       \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    except Exception as e:
      raise GrepError('Encountered error while grepping for lattice vectors '  \
                      'in file {1.directory}/{1.STDOUT}: '                     \
                      '{0.__class__.__name__} -- {0}'.format(e, self) )
    finally: file.close()

  @property
  @make_cached
  def symmetry_operators(self):
    """ Greps symmetry operators from file. """
    from numpy import array, zeros, dot
    try:
      file = self.__stdout__()
      for line in file:
        if line.split() == ['V', 'INV', 'ROTATION', 'MATRICES', 'TRANSLATOR']:
          break
      symops = []
      for line in file:
        data = line.split()
        if len(data) != 14: break
        symops.append(zeros((4,3), dtype='float64'))
        symops[-1][:3] = array(data[2:11], dtype='float64').reshape(3,3).T
        symops[-1][3] = array(data[-3:], dtype='float64')
        symops[-1][3] = dot(self.structure.cell, symops[-1][3])
      return array(symops)
    except Exception as e:
      raise GrepError('Encountered error while grepping for sym ops '          \
                      'in file {1.directory}/{1.STDOUT}: '                     \
                      '{0.__class__.__name__} -- {0}'.format(e, self) )
    finally: file.close()

  @property
  @make_cached
  def functional(self):
    """ Reads functional from output file. """
    from .functional import Functional
    from .parse import parse
    with self.__stdout__() as file: b = parse(file)
    result = Functional()
    result.read_input(b)
    return result

  @property
  @make_cached
  def title(self):
    """ Title of the calculations. """
    search = " EEEEEEEEEE STARTING  DATE"
    with self.__stdout__() as file: 
      for line in file: 
        if len(line) <= len(search): continue
        if line[:len(search)] == search: break
      try: line = file.next().rstrip().lstrip()
      except StopIteration: return ""
      else: return line
      
  @property
  @make_cached
  def start_date(self):
    """ Title of the calculations. """
    from datetime import datetime
    search = " EEEEEEEEEE STARTING  DATE"
    with self.__stdout__() as file: 
      for line in file: 
        if len(line) <= len(search): continue
        if line[:len(search)] == search: break
      line = line.split()[3:]
      day, month, year = [int(u) for u in line[:3]]
      time = line[4].split(':')
      hour, minute, second = int(time[0]), int(time[1]), time[2]
      second = int(second[:second.find('.')])
      return datetime( year=year, month=month, day=day,
                       hour=hour, minute=minute, second=second )

  @property
  @make_cached
  def end_date(self):
    """ Title of the calculations. """
    from datetime import datetime
    from ..error import GrepError
    pattern = "E+\s+TERMINATION\s+DATE\s*(\d+)\s+(\d+)\s+(\d+)\s+TIME\s+(\S+)"
    regex = self._find_last_STDOUT(pattern)
    if regex is None: 
      raise GrepError( 'Could not grep end time from '                         \
                       '{0.directory}/{0.STDOUT}.'.format(self) )

    day    = int(regex.group(1))
    month  = int(regex.group(2))
    year   = int(regex.group(3))
    hour, minute, second = regex.group(4).split(':')
    hour   = int(hour)
    minute = int(minute)
    second = int(second[:second.find('.')])
    return datetime( year=year, month=month, day=day,
                     hour=hour, minute=minute, second=second )

  @property
  @make_cached
  def cputime(self):
    """ Total CPU time. """
    from datetime import timedelta
    from ..error import GrepError
    regex = self._find_first_STDOUT("TOTAL\s+CPU\s+TIME\s+=\s+(\S+)")
    if regex is None: raise GrepError("Could not grep total cpu time.")
    return timedelta(seconds=float(regex.group(1)))

  @property
  @make_cached
  def moqgad_time(self):
    """ MOQGAD time. """
    from datetime import timedelta
    from re import compile
    regex = compile('T+\s+MOQGAD\s+TELAPSE\s+(\S+)\s+TCPU\s+(\S+)')
    with self.__stdout__() as file:
      results = []
      for u in regex.finditer(file.read()): 
        results.append(timedelta(seconds=float(u.group(2))))
      return results

  @property
  @make_cached
  def shellx_time(self):
    """ SHELLX time. """
    from datetime import timedelta
    from re import compile
    regex = compile('T+\s+SHELLX\s+TELAPSE\s+(\S+)\s+TCPU\s+(\S+)')
    with self.__stdout__() as file:
      results = []
      for u in regex.finditer(file.read()): 
        results.append(timedelta(seconds=float(u.group(2))))
      return results

  @property
  @make_cached
  def monmo3_time(self):
    """ MONMO3 time. """
    from datetime import timedelta
    from re import compile
    regex = compile('T+\s+MONMO3\s+TELAPSE\s+(\S+)\s+TCPU\s+(\S+)')
    with self.__stdout__() as file:
      results = []
      for u in regex.finditer(file.read()): 
        results.append(timedelta(seconds=float(u.group(2))))
      return results

  @property
  @make_cached
  def numdft_time(self):
    """ NUMDFT time. """
    from datetime import timedelta
    from re import compile
    regex = compile('T+\s+NUMDFT\s+TELAPSE\s+(\S+)\s+TCPU\s+(\S+)')
    with self.__stdout__() as file:
      results = []
      for u in regex.finditer(file.read()): 
        results.append(timedelta(seconds=float(u.group(2))))
      return results

  @property
  @make_cached
  def fdik_time(self):
    """ FDIK time. """
    from datetime import timedelta
    from re import compile
    regex = compile('T+\s+FDIK\s+TELAPSE\s+(\S+)\s+TCPU\s+(\S+)')
    with self.__stdout__() as file:
      results = []
      for u in regex.finditer(file.read()): 
        results.append(timedelta(seconds=float(u.group(2))))
      return results

  @property
  @make_cached
  def pdig_time(self):
    """ PDIG time. """
    from datetime import timedelta
    from re import compile
    regex = compile('T+\s+PDIG\s+TELAPSE\s+(\S+)\s+TCPU\s+(\S+)')
    with self.__stdout__() as file:
      results = []
      for u in regex.finditer(file.read()): 
        results.append(timedelta(seconds=float(u.group(2))))
      return results

  @property
  @make_cached
  def input_title(self):
    """ Title as given on input. """
    from .parse import parse
    with self.__stdout__() as file: tree = parse(file)
    return tree.keys()[0]

  @property
  @make_cached
  def is_molecular(self):
    """ True if a molecular calculation """
    pattern = "^\s+(MOLECULAR|CRYSTAL)\s+CALCULATION\s*$"
    regex = self._find_last_STDOUT(pattern)
    if regex is None:
      raise GrepError('Could not determine whether molecular calculation')
    return regex.group(1) == 'MOLECULAR'
    

  def _parsed_tree(self):
    """ Returns parsed input tree. """
    from .parse import parse
    with self.__stdout__() as file: tree = parse(file)
    title = tree.keys()[0]
    return tree[title]

  @property
  @make_cached
  def input_crystal(self):
    """ Input structure, CRYSTAL format. 
    
        CRYSTAL_ defines structure in a functional way, starting from an input
        structure on which are applied a sequence of transformations. This
        structure format is grepped directly from the output.
    """
    from .crystal import Crystal
    from .molecule import Molecule
    from .. import CRYSTAL_geom_blocks as starters
    from ..error import IOError, NotImplementedError

    tree = self._parsed_tree()
    found = False
    for starter in starters:
      if starter in tree.keys(): found = True; break
    if found == False:
      raise IOError('Could not find start of input in file.')
    if starter.lower() != 'crystal': 
      raise NotImplementedError('Can only read 3d structures.')

    result = Molecule() if self.is_molecular else Crystal()
    result.read_input(tree[starter])
    return result

  def _grep_structure(self, file):
    """ Greps structure from file.

        This method finds the next structure from the file. It can be used in
        conjunction with some search method and seek to extract either the
        first or the last structure, or anything in between.
    """
    from numpy import array, identity
    from ..crystal import Structure
    from ..error import GrepError
    from .basis import specie_name
    result = Structure()
    try: 
      file.next(); file.next() # move to first line.
      for line in file:
        line = line.split()
        if len(line) != 7: break
        type = specie_name(int(line[2]))
        asymmetric = line[1] == 'T' 
        result.add_atom( pos=array(line[4:7], dtype='float64'),
                         type=type, label=int(line[0]), asymmetric=asymmetric,
                         group=line[3] )
        
      # If a molecule, set cell to 500.0 as in CRYSTAL
      if self.is_molecular:
        self.cell = identity(3, dtype='float64') * 500.0
      else:
        # then find cell.
        header = 'DIRECT LATTICE VECTORS CARTESIAN '                             \
                 'COMPONENTS (ANGSTROM)'.split()
        for line in file: 
          line = line.split()
          if len(line) != 6: continue
          if line == header: break
        try: file.next()
        except StopIteration: raise GrepError('File is incomplete.')
        result.cell = array( [file.next().split() for i in xrange(3)],
                             dtype='float64' )
      
        # Then re-reads atoms, but in cartesian coordinates.
        for i in xrange(6): file.next()
        for atom in result:
          atom.pos = array(file.next().split()[3:6], dtype='float64')
    except StopIteration: raise GrepError('Unexpected end of file')

    # adds more stuff
    try: title = self.title
    except: pass
    else:
     if len(title) > 0: result.name = title
    return result

  def _find_structure(self, file):
    """ Yields positions in file where structure starts. """
    # first finds atoms -- the one with whether they are in the asymmetric
    # unit. The file is scanned to the end unit so the last position of the
    # header can be found. The, we use seek to go back to that position  and
    # find the result. We cannot do this with file.next() because of
    # read-ahead buffering issues. Also, we go back one line so things work better.
    header = 'ATOMS IN THE ASYMMETRIC UNIT'.split()
    while True:
      line = file.readline()
      if not line: break
      line = line.split()
      if len(line) < 5: continue
      if line[:5] == header: yield file.tell()

  @property
  @make_cached
  def input_structure(self):
    """ Input structure, LaDa format. """
    with self.__stdout__() as file:
      pos = None
      for pos in self._find_structure(file): break
      if pos is None: raise GrepError('Could extract structure from file')
      file.seek(pos, 0)
      return self._grep_structure(file)


  @property
  @make_cached
  def structure(self):
    """ Output structure, LaDa format. """
    with self.__stdout__() as file:
      pos = None
      for pos in self._find_structure(file): continue
      if pos is None: raise GrepError('Could not extract structure from file')
      file.seek(pos, 0)
      return self._grep_structure(file)

  @property
  @make_cached
  def crystal(self):
    """ Output structure, CRYSTAL format. 
    
        Creates a structure from the output, including any atomic and
        cell-shape movements.
    """
    from numpy import dot, identity, abs, array, any, all
    from numpy.linalg import inv, det
    from ..error import internal 
    from .geometry import DisplaceAtoms, Elastic
    
    # Check whether this is a geometry optimization run.
    if self._find_first_STDOUT("STARTING GEOMETRY OPTIMIZATION") is None:
      return self.input_crystal

    # Find last operation which is neither ELASTIC nor ATOMDISP
    incrys, instruct, i = self.input_crystal, self.input_structure, 0
    looped = False
    for i, op in enumerate(incrys[::-1]):
      if op.keyword.lower() not in ['elastic', 'atomdisp']: break
      looped = True

    # deduce structure - last changes in cell-shape or atomic displacements.
    if looped:
      if incrys[-i].keyword.lower() in ['elastic', 'atomdisp']: i += 1
      incrys = incrys.copy()
      incrys[:] = incrys[:len(incrys)-i]
      instruct = incrys.eval()

    # create symmetric strain
    inv_in = inv(instruct.cell)
    epsilon = dot(self.structure.cell, inv_in) - identity(3, dtype='float64')
    epsilon = 0.5 * (epsilon + epsilon.T)
    cell = dot(identity(3) + epsilon, instruct.cell)
    inv_out = inv(cell)
    if any(abs(cell - self.structure.cell) > 1e-8):
      raise internal('Could not create symmetric strain matrix')

    # create field displacement
    field = [ dot(cell, dot(inv_out, a.pos) - dot(inv_in, b.pos))
              for a, b in zip(self.structure, instruct)
              if a.asymmetric ]
    field = array(field)

    # check if changes:
    if all(abs(epsilon) < 1e-8) and all(abs(field.flatten()) < 1e-8): 
      return self.input_crystal

    result = incrys.copy()
    # add cell shape changes
    if any(abs(epsilon) > 1e-8): 
      a = Elastic()
      a.matrix = epsilon
      a.is_epsilon = True
      a.const_volume = abs(det(cell) - det(instruct.cell)) < 1e-8
      result.append(a)
    # Add displacements 
    if any(abs(field.flatten()) > 1e-8):
      a = DisplaceAtoms(keepsymm=True)
      atoms = [u for u in instruct if u.asymmetric]
      for atom, disp in zip(atoms, field):
        if any(abs(disp) > 1e-8): a.add_atom(type=atom.label, pos=disp)
      result.append(a)
    return result


  @property
  @make_cached
  def total_energy(self):
    """ Total energy. """
    from quantities import hartree
    pattern = "TOTAL ENERGY\(\S+\)\(AU\)\(\s*\d+\)\s*(\S+)\s*DE"
    regex = self._find_last_STDOUT(pattern)
    if regex is None: 
      raise GrepError( 'Could not grep total energy from '                     \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    return float(regex.group(1)) * hartree

  @property
  @make_cached
  def total_energies(self):
    """ Total energies until convergences. """
    from numpy import array
    from quantities import hartree
    pattern = r"CYC\s+(?:\d+)\s+ETOT\(AU\)\s+(\S+)"
    result = [u.group(1) for u in self._search_STDOUT(pattern)]
    return array(result, dtype='float64') * hartree

  @property
  @make_cached
  def delta_energy(self):
    """ Difference in energy between last two minimization steps. """
    from quantities import hartree
    pattern = "TOTAL ENERGY\(\S+\)\(AU\)\(\s*\d+\)\s*\S+\s*DE(\S+)\s*tester"
    regex = self._find_last_STDOUT(pattern)
    if regex is None: 
      raise GrepError( 'Could not grep delta energy from '                     \
                       '{0.directory}/{0.STDOUT}.'.format(self) )
    return float(regex.group(1)) * hartree

  def iterfiles(self, **kwargs):
    """ iterates over input/output files. 
    
        :param bool input: Include INCAR file
        :param bool wavefunctions: Include WAVECAR file
        :param bool structure: Include POSCAR file
    """
    from os.path import exists, join
    files = [self.STDOUT]
    if kwargs.get('input', False):   files.append('crystal.d12')
    if kwargs.get('wavefunctions', False): files.append('crystal.f98')
    if kwargs.get('structure', False):  files.append('crystal.34')
    for file in files:
      file = join(self.directory, file)
      if exists(file): yield file

  @property
  @make_cached
  def optgeom_convergence(self): 
    """ True if optgeom convergence was achieved

        False if it was not achieved, and None in all other cases.
        Also returns None if this was not a geometry optimization run.
    """
    try: 
      if self._find_first_STDOUT('CONVERGENCE TESTS SATISFIED') is not None:
        return True
      if self._find_first_STDOUT('CONVERGENCE ON GRADIENTS SATISFIED') is not None:
        return True
      if self._find_first_STDOUT('CONVERGENCE TESTS UNSATISFIED') is not None:
        return False
    except: pass
    return None

  @property
  @make_cached
  def surface_area(self):
    """ Surface area of 2d slabs. """
    from quantities import angstrom
    from ..error import GrepError
    result = self._find_last_STDOUT('AREA\s+OF\s+THE\s+2D\s+CELL\s+(\S+)')
    if result is None: raise GrepError('Likely not a 2D slab')
    return float(result.group(1)) * angstrom * angstrom 

  @property
  @make_cached
  def slab_height(self):
    """ Height of the slab. """
    from quantities import angstrom
    from ..error import GrepError
    a = self.surface_area.magnitude
    result = self._find_last_STDOUT('VOLUME\s+OF\s+THE\s+3D\s+CELL\s+(\S+)')
    if result is None: raise GrepError('Likely not a 2D slab')
    return float(result.group(1)) / a * angstrom 


  @property
  @make_cached
  def optgeom_iterations(self):
    """ Number of geometric iterations. """
    a = self.optgeom_convergence
    if a is None: return 0
    result = self._find_last_STDOUT('POINTS\s+(\d+)\s+\*')
    return int(result.group(1))

class Extract(AbstractExtractBase, OutputSearchMixin, ExtractBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, **kwargs):
    """ Initializes extraction object. 
    
        :param directory: 
          Directory where the OUTCAR resides. 
          It may also be the path to an OUTCAR itself, if the file is not
          actually called OUTCAR.
    """
    from os.path import exists, isdir, basename, dirname
    from lada.misc import RelativePath
       
    self.stdout = 'crystal.out'
    """ Name of file to grep. """
    if directory is not None:
      directory = RelativePath(directory).path
      if exists(directory) and not isdir(directory):
        self.STDOUT = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractBase.__init__(self)
    OutputSearchMixin.__init__(self)

  @property
  def success(self):
    try: self.end_date
    except: 
      return self.optgeom_convergence is not None
    return True

class MassExtract(AbstractMassExtract):
  """ Extracts all CRYSTAL calculations in directory and sub-directories. 
    
      Trolls through all subdirectories for vasp calculations, and organises
      results as a dictionary where keys are the name of the diretories.

      Usage is simply:

      >>> from lada.dftcrystal import MassExtract
      >>> a = MassExtract('path') # or nothing if path is current directory.
      >>> a.success
      {
        '/some/path/':      True,
        '/some/other/path': True
      }
  """
  DefaultExtract = Extract
  """ Extraction class to use by default. 

      This is mostly to allow for some kind of monkey patching.
  """
  def __init__(self, path = None, glob=None, **kwargs):
    """ Initializes MassExtract.
    
    
        :param str path:
            Root directory for which to investigate all subdirectories.
            If None, uses current working directory.
        :param str glob:
            If not None, then uses this glob to find the files to grep. These
            should be CRYSTAL output files.
        :param kwargs:
            Passed on to
            :py:class:`~lada.jobfolder.extract.AbstractMassExtract`
    """
    from os import getcwd
    if path is None: path = getcwd()
    # this will throw on unknown kwargs arguments.
    super(MassExtract, self).__init__(path=path, **kwargs)
    self.glob = glob
    """ Pattern from which to determine output file. """

  def __iter_alljobs__(self):
    """ Goes through all directories with an OUTVAR. """
    from os import walk
    from os.path import relpath, join
    from glob import iglob
    from .relax import RelaxExtract

    if self.glob is None:
      for dirpath, dirnames, filenames in walk( self.rootpath, topdown=True, 
                                                followlinks=True ):
        if 'crystal.out' not in filenames: continue
        if 'relax' in dirnames:
          dirnames[:] = [u for u in dirnames if u != 'relax']
          try: result = RelaxExtract(join(self.rootpath, dirpath))
          except:
            try: result = self.__class__.DefaultExtract(join(self.rootpath, dirpath))
            except: continue
        else: 
          try: result = self.__class__.DefaultExtract(join(self.rootpath, dirpath))
          except: continue
        
        yield join('/', relpath(dirpath, self.rootpath)), result
    else: 
      for dirpath, dirnames, filenames in walk( self.rootpath, topdown=True, 
                                                followlinks=True ):
        for filename in iglob(join(join(self.rootpath, dirpath), self.glob)):
          try: result = self.__class__.DefaultExtract(filename)
          except: continue
          else: yield join('/',relpath(filename, self.rootpath)), result

  def __copy__(self):
    """ Returns a shallow copy. """
    result = self.__class__(self.rootpath)
    result.__dict__.update(self.__dict__)
    return result

  @property
  def _attributes(self): 
    """ Returns __dir__ set special to the extraction itself. """
    return list(set( [ u for u in dir(self.DefaultExtract) if u[0] != '_' ] 
                     + ['details'] ))
