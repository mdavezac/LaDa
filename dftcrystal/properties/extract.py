""" Subpackage containing extraction methods for CRYSTAL's property output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract', 'MassExtract']

from ...tools import make_cached
from ...tools.extract import search_factory, AbstractExtractBase
from ...error import GrepError
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'prop.out')

class ExtractBase(object):
  """ Implementation class for extracting CRYSTAL's property output. """
  def __init__(self):
    """ Initializes the extraction class. """
    super(ExtractBase, self).__init__()

  @property
  @make_cached
  def is_bandstructure(self):
    """ True if contains band-structure calculation. """
    result = self._find_first_STDOUT("""^\s*\*\s*BAND\s+STRUCTURE\s*\*\s*$""")
    return result is not None

  @property
  @make_cached
  def end_date(self):
    """ Title of the calculations. """
    from datetime import datetime
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
  def bandstructure(self):
    """ Computed band-structure. 

        If the band-structure was *not* computed, returns None.
        If it was, this is a named tuple ``(kpoints=[...], eigenvalues=[...])``.
        The kpoints are in the same order as per the actual calculation. They
        are always in fractional coordinates. None of that fake floating point
        iss thingie from crystal.  The eigenvalues are given in hartree, as a
        ``kpoint x band`` matrix.
    """
    from collections import namedtuple
    from re import M
    from numpy import array, sum, sqrt
    # if not a band-structure, then returns.
    if not self.is_bandstructure: return None

    # create the result class.
    BandStructure = namedtuple('BandStructure', ['kpoints', 'eigenvalues'])

    # then greps k-points.
    regex = r"\s*CARTESIAN\s+COORD\.\s+\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)"       \
            r"\s*:\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)\s*STEP\s+(\S+)\s*$"
    kpoints = []
    for kpre in self._search_STDOUT(regex, M):
      start = array([kpre.group(i) for i in xrange(1, 4)], dtype='float64')
      end = array([kpre.group(i) for i in xrange(4, 7)], dtype='float64')
      step = float(kpre.group(7))
      print start, end, step, sqrt(sum((start-end)**2)) / step



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
       
    self.STDOUT = 'prop.out'
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
    if not exists(join(self.directory, self.STDOUT)): return False
    try: self.end_date
    except: return False
    else: return True

