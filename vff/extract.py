""" Subpackage containing extraction methods for CRYSTAL output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']

from ..tools import make_cached
from ..tools.extract import search_factory, AbstractExtractBase
from ..error import GrepError
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'vff.out')

class ExtractBase(object):
  """ Implementation class for extracting data from VFF output """

  def __init__(self):
    """ Initializes the extraction class. """
    super(ExtractBase, self).__init__()

  @property
  @make_cached
  def start_date(self):
    """ Starting date of the calculation. """
    from datetime import datetime
    regex = self._find_first_STDOUT('^Start date:\s*(.*)$')
    if regex is None: raise GrepError('Could not find starting date.')
    regex = regex.group(1).split()
    regex = list(regex[0].split('-')) + list(regex[1].split(':'))
    regex = regex[:-1] + list(regex[-1].split('.'))
    regex = [int(u) for u in regex]
    return datetime(*regex)
  @property
  @make_cached
  def end_date(self):
    """ Starting date of the calculation. """
    from datetime import datetime
    regex = self._find_first_STDOUT('^End date:\s*(.*)$')
    if regex is None: raise GrepError('Could not find end date.')
    regex = regex.group(1).split()
    regex = list(regex[0].split('-')) + list(regex[1].split(':'))
    regex = regex[:-1] + list(regex[-1].split('.'))
    regex = [int(u) for u in regex]
    return datetime(*regex)



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
       
    self.STDOUT = 'vff.out'
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
    from os.path import exists, join
    try: self.end_date
    except: return False
    return True

  def iterfiles(self, **kwargs):
    """ iterates over input/output files. 
    
        :param bool input: Include INCAR file
        :param bool wavefunctions: Include WAVECAR file
        :param bool structure: Include POSCAR file
    """
    from os.path import exists, join
    from .. import CRYSTAL_filenames as filenames
    file = join(self.directory, self.STDOUT)
    if exists(file): yield file

del make_cached
del search_factory
