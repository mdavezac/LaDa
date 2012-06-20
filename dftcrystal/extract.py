from lada.functools import make_cached
from lada.functools.extract import search_factory, AbstractExtractBase
from lada.error import GrepError
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__)

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
        self.stdout = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractBase.__init__(self)
    OutputSearchMixin.__init__(self)

  @property
  def success(self): return True
