""" Subpackage containing extraction methods for vasp parameters from vasp output. 

    Extaction objects are implemented as a mix and mash of bases classes. The
    reason for this is we want to isolate functionality specific to DFT and GW,
    and specific to reading *real* OUTCAR files and *database* OUTCAR files. 
"""
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from ...functools.extract import AbstractExtractBase
from ._common import Extract as ExtractCommonBase
from ._dft import Extract as ExtractDFTBase
from ._mixin import IOMixin

class ExtractCommon(AbstractExtractBase, ExtractCommonBase, IOMixin):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, **kwargs):
    """ Initializes extraction object.
    
    
        :Parameters:
          directory : str or None
            Path to OUTCAR file. Can also be the directory if the OUTCAR is
            named "OUTCAR".
    """
    from os.path import exists, isdir, basename, dirname
    from ...misc import RelativePath
       
    outcar = None
    if directory is not None:
      directory = RelativePath(directory).path
      if exists(directory) and not isdir(directory):
        outcar = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractCommonBase.__init__(self)
    IOMixin.__init__(self, directory, OUTCAR=outcar, **kwargs)

  @property
  def success(self):
    """ True if calculation was successfull. """
    return ExtractCommonBase.success.__get__(self)

class ExtractDFT(AbstractExtractBase, IOMixin, ExtractDFTBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, **kwargs):
    """ Initializes extraction object. """
    from os.path import exists, isdir, basename, dirname
    from ...misc import RelativePath
       
    outcar = None
    if directory is not None:
      directory = RelativePath(directory).path
      if exists(directory) and not isdir(directory):
        outcar = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractDFTBase.__init__(self)
    IOMixin.__init__(self, directory, OUTCAR=outcar, **kwargs)
  @property
  def success(self):
    """ True if calculation was successfull. """
    return ExtractDFTBase.success.__get__(self)

def Extract(directory=None, **kwargs):
  """ Chooses between DFT or GW extraction object, depending on OUTCAR.
  
        :Parameters:
          outcar : str or None
            Path to OUTCAR file. Can also be the directory if the OUTCAR is
            named "OUTCAR". Defaults to None, in which case it uses the current
            working directory.
  """
  from os import getcwd
  from os.path import exists, isdir, join
  from ...misc import RelativePath
  # checks for GW calculations.
  directory = getcwd() if directory is None else RelativePath(directory).path
  if exists(join(directory, 'GW')):
    from os.path import basename
    from glob import iglob
    from re import match
    outcar = join(directory, 'GW')
    gwiters = [ int(basename(u)) for u in iglob(join(outcar, "[0-9]*"))\
                if match(r"\d+", basename(u)) and isdir(u) ]
    if len(gwiters) > 0: outcar = join(outcar, str(max(gwiters)))

  result = ExtractCommon(directory=directory, **kwargs)
  if result.is_dft: return ExtractDFT(directory=directory, **kwargs) 
# elif result.is_gw: return ExtractGW(directory=directory, **kwargs) 
  return result

