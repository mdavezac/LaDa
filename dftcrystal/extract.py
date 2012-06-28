from lada.functools import make_cached
from lada.functools.extract import search_factory, AbstractExtractBase
from lada.error import GrepError
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
    if result is None: raise GrepError('Could not grep crystal family')
    return result.group(1).lower()
  @property 
  @make_cached
  def crystal_class(self):
    """ Crystal class """
    result = self._find_first_STDOUT(r"^\s*CRYSTAL CLASS\s*(?:.+):\s*(.+)$")
    if result is None: raise GrepError('Could not grep crystal class')
    return result.group(1).lower().rstrip().lstrip()

  @property 
  @make_cached
  def centrosymmetric(self):
    """ Crystal class """
    result = self._find_first_STDOUT(r"^\s*SPACE GROUP\s*\((\S+)\)\s*:")
    if result is None: raise GrepError('Could not grep centro-symmetricity')
    return result.group(1).lower().rstrip().lstrip() != 'noncentrosymmetric'

  @property 
  @make_cached
  def space_group(self):
    """ Crystal class """
    result = self._find_first_STDOUT(r"^\s*SPACE GROUP\s*(?:.+)\s*:\s*(.+)")
    if result is None: raise GrepError('Could not grep space-group')
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
      raise GrepError('Reached end of file whend grepping '                    \
                      'direct lattice cell vectors.')
    except Exception as e:
      raise GrepError('Encountered error while grepping for lattice vectors: ' \
                      '{0.__class__.__name__} -- {0}'.format(e) )
    finally: file.close()

  @property
  @make_cached
  def symmetry_operators(self):
    """ Greps symmetry operators from file. """
    from numpy import array, zeros
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
        symops[-1][3, :] = array(data[-3:], dtype='float64')
      return array(symops)
    except Exception as e:
      raise GrepError('Encountered error while grepping for sym ops:'          \
                      '{0.__class__.__name__} -- {0}'.format(e) )
    finally: file.close()

  @property
  @make_cached
  def functional(self):
    """ Reads functional from output file. """
    from .functional import read_input
    with self.__stdout__() as file: return read_input(file)

  @property
  @make_cached
  def title(self):
    """ Title of the calculations. """
    search = " EEEEEEEEEE STARTING  DATE"
    with self.__stdout__() as file: 
      for line in file: 
        if len(line) <= len(search): continue
        if line[:len(search)] == search: break
      line = file.next().rstrip().lstrip()
      if line[-1] == '\n': line = line[:-1]
      return line
      
  @property
  @make_cached
  def date(self):
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
  def input_title(self):
    """ Title as given on input. """
    from .parse import parse
    with self.__stdout__() as file: tree = parse(file)
    return tree.keys()[0]

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
    from .geometry import Crystal
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

    result = Crystal()
    result.read_input(tree[starter])
    return result

  def _grep_structure(self, file):
    """ Greps structure from file.

        This method finds the next structure from the file. It can be used in
        conjunction with some search method and seek to extract either the
        first or the last structure, or anything in between.
    """
    from numpy import array
    from ..crystal import Structure
    result = Structure()
    file.next(); file.next() # move to first line.
    for line in file:
      line = line.split()
      if len(line) != 7: break
      type = int(line[2])
      if type < 100: type = line[3]
      asymmetric = line[1] == 'T' 
      result.add_atom( pos=array(line[4:7], dtype='float64'),                 \
                       type=type, label=int(line[0]), asymmetric=asymmetric )
      
    # then find cell.
    header = 'DIRECT LATTICE VECTORS CARTESIAN '                             \
             'COMPONENTS (ANGSTROM)'.split()
    for line in file: 
      line = line.split()
      if len(line) != 6: continue
      if line == header: break
    file.next()
    result.cell = array( [file.next().split() for i in xrange(3)],           \
                         dtype='float64' )

    # Then re-reads atoms, but in cartesian coordinates.
    for i in xrange(6): file.next()
    for atom in result:
      atom.pos = array(file.next().split()[3:6], dtype='float64')

    # adds more stuff
    result.name = self.title
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
      pos = self._find_structure(file).next()
      file.seek(pos, 0)
      return self._grep_structure(file)


  @property
  @make_cached
  def structure(self):
    """ Input structure, LaDa format. """
    with self.__stdout__() as file:
      for pos in self._find_structure(file): continue
      file.seek(pos, 0)
      return self._grep_structure(file)

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
       
    self.stdout = 'stdout'
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
  def success(self): return True
