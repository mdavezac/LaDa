""" Subpackage containing extraction methods for vasp parameters from vasp output. 

    Extaction objects are implemented as a mix and mash of bases classes. The
    reason for this is we want to isolate functionality specific to DFT and GW,
    and specific to reading *real* OUTCAR files and *database* OUTCAR files. 
"""
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from ...opt import AbstractExtractBase, OutcarSearchMixin as SearchMixin
from ._common import Extract as ExtractCommonBase
from ._dft import Extract as ExtractDFTBase
from ._gw import Extract as ExtractGWBase
from ._mixin import IOMixin

class ExtractCommon(AbstractExtractBase, ExtractCommonBase, IOMixin, SearchMixin):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, comm=None, **kwargs):
    """ Initializes extraction object. """
    from os.path import exists, isdir, basename, dirname
    # checks if path or directory
    if directory != None and exists(directory) and not isdir(directory):
      kwargs['OUTCAR'] = basename(directory)
      directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory, comm)
    ExtractCommonBase.__init__(self)
    IOMixin.__init__(self, directory, **kwargs)
    SearchMixin.__init__(self)

  @property
  def success(self):
    """ True if calculation was successfull. """
    return ExtractCommonBase.success.__get__(self)

class ExtractDFT(ExtractCommon, ExtractDFTBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, directory, comm, **kwargs)
    ExtractDFTBase.__init__(self)

class ExtractGW(ExtractCommon, ExtractGWBase):
  """ Extracts GW data from an OUTCAR. """
  def __init__(self, directory=None, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, directory, comm, **kwargs)
    ExtractGWBase.__init__(self)

def Extract(*args, **kwargs): 
  """ Chooses between DFT or GW extraction object, depending on OUTCAR. """
  a = ExtractCommon(*args, **kwargs)
  try: which = ExtractDFT if a.is_dft else ExtractGW
  except: which = ExtractCommon
  return which(*args, **kwargs)

def ExtractGW_deprecated(*args, **kwargs):
  """ Deprecated. Please use vasp.Extract instead. """
  from warnings import warn
  warn('ExtractGW is deprecated. Please use vasp.Extract instead.', DeprecationWarning)
  return Extract(*args, **kwargs)
    
try: from ... import jobs
except ImportError: pass
else: 
  __all__.append('MassExtract')
  class MassExtract(jobs.AbstractMassExtract):
    """ Propagates vasp extractors from all subdirectories.
    
        Trolls through all subdirectories for vasp calculations, and organises
        results as a dictionary where keys are the name of the diretory.
    """
    def __init__(self, path = None, Extract = None, **kwargs):
      """ Initializes MassExtract.
      
      
          :Parameters:
            path : str or None
              Root directory for which to investigate all subdirectories.
              If None, uses current working directory.
            Extract : `lada.vasp.Extract`
              Extraction class to use. 
            kwargs : dict
              Keyword parameters passed on to AbstractMassExtract.

          :kwarg naked_end: True if should return value rather than dict when only one item.
          :kwarg unix_re: converts regex patterns from unix-like expression.
      """
      from os import getcwd
      from os.path import exists, isdir
      from . import Extract as VaspExtract
      from ...opt import RelativeDirectory

      # this will throw on unknown kwargs arguments.
      jobs.AbstractMassExtract.__init__(self, **kwargs)

      self.Extract = Extract if Extract != None else VaspExtract
      """ Extraction class to use. """

      if path == None: path = getcwd()
      self._rootdir = RelativeDirectory(path, hook=self.uncache)
      """ Root of the directory-tree to trawl for OUTCARs. """
      
      self.OUTCAR = "OUTCAR"
      """ Name of the OUTCAR file. """
      assert exists(self.rootdir), RuntimeError("Path {0} does not exist.".format(self.rootdir))
      assert isdir(self.rootdir), RuntimeError("Path {0} is not a directory.".format(self.rootdir))

    @property
    def rootdir(self): 
      """ Root of the directory-tree to trawl for OUTCARs. """
      return self._rootdir.path
    @rootdir.setter
    def rootdir(self, value): self._rootdir.path = value

    def __iter_alljobs__(self):
      """ Goes through all directories with an OUTVAR. """
      from os import walk, getcwd
      from os.path import abspath, relpath, abspath, join

      for dirpath, dirnames, filenames in walk(self.rootdir, topdown=True, followlinks=True):
        if self.OUTCAR not in filenames: continue

        try: result = self.Extract(join(self.rootdir, dirpath))
        except: continue

        result.OUTCAR = self.OUTCAR
        yield join('/', relpath(dirpath, self.rootdir)), result

    def __copy__(self):
      """ Returns a shallow copy. """
      result = self.__class__(self.rootdir)
      result.__dict__.update(self.__dict__)
      return result

      


