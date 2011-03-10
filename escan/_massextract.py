""" Module to extract esca and vff ouput. """
__docformat__ = "restructuredtext en"
__all__ = ['MassExtract']

from ..jobs import AbstractMassExtractDirectories
from ._extract import Extract


class MassExtract(AbstractMassExtractDirectories):
  """ Extracts all escan calculations nested within a given input directory. """
  def __init__(self, path = ".", **kwargs):
    """ Initializes AbstractMassExtractDirectories.
    
    
        :Parameters:
          path : str or None
            Root directory for which to investigate all subdirectories.
            If None, uses current working directory.
          kwargs : dict
            Keyword parameters passed on to AbstractMassExtractDirectories.

        :kwarg naked_end: True if should return value rather than dict when only one item.
        :kwarg unix_re: converts regex patterns from unix-like expression.
    """
    # this will throw on unknown kwargs arguments.
    if 'Extract' not in kwargs: kwargs['Extract'] = Extract
    super(MassExtract, self).__init__(path, **kwargs)
    del self.__dict__['Extract']

  def __iter_alljobs__(self):
    """ Goes through all directories with a contcar. """
    from os import walk, getcwd
    from os.path import relpath, join, dirname
    from ._extract import Extract as EscanExtract

    OUTCAR = EscanExtract().OUTCAR
    self.__dict__['__attributes'] = set()
    for dirpath, dirnames, filenames in walk(self.rootdir, topdown=True, followlinks=True):
      i = len(dirnames) - 1
      for dirname in list(dirnames[::-1]):
        if dirname in ['AE', 'CBM', 'VBM']: dirnames.pop(i)
        elif len(dirname) > 7 and dirname[:7] == 'kpoint_': dirnames.pop(i)
        i -= 1

      try: result = self.Extract(join(self.rootdir, dirpath), comm = self.comm)
      except: continue

      if result.__class__.__name__ != 'ExtractEscanFail': 
        yield join('/', relpath(dirpath, self.rootdir)), result
        self.__dict__['__attributes'] |= set([u for u in dir(result) if u[0] != '_'])
    self.__dict__['__attributes'] = list(self.__dict__['__attributes'])
  @property
  def _attributes(self):
    """ List of attributes inherited from escan extractors. """
    if '__attributes' not in self.__dict__:
      self.__class__._extractors.__get__(self)
    return self.__dict__['__attributes']

  def __is_calc_dir__(self, dirpath, dirnames, filenames):
    """ Returns true this directory contains a calculation. """
    from os.path import dirname
    from ._extract import Extract as EscanExtract
    dir = dirpath[len(dirname(dirpath)):]
    if dir in ['AE', 'VBM', 'CBM']: return False
    if len(dir) > 7 and dir[:7] == 'kpoint_': return False
    OUTCAR = EscanExtract().OUTCAR
    return OUTCAR in filenames

  @staticmethod
  def Extract(directory=None, **kwargs):
    """ Extracts from directory for Escan, Kescan, and bandgap. 
    
        Looks at file and directory structure to deduce whether these are simple
        escan calculations, kescan calculations, or bandgap calculations.
    """
    from os import getcwd
    from os.path import exists, join
    from ..opt import RelativeDirectory
    from ._extract import Extract
    from .kescan import Extract as KExtract
    from ._bandgap import extract as extract_bg
    OUTCAR = kwargs.get('OUTCAR', Extract().OUTCAR)
    comm = kwargs.get('comm', None)

    class ExtractEscanFail(object):
      @property
      def success(self): return False
    if directory == None: directory = getcwd()
    else: directory = RelativeDirectory(directory).path
    if not exists(directory): return ExtractEscanFail()

    # found bandgap calculation.
    if exists(join(directory, 'AE')): return extract_bg(directory, **kwargs)
    if not exists(join(directory, OUTCAR)): return ExtractEscanFail()
    if exists(join(directory, 'CBM')) and exists(join(directory, 'VBM')):
      return extract_bg(directory, **kwargs)

    # found escan or kescan calculation.
    result = Extract(directory, **kwargs)
    return result if result.functional.__class__.__name__ != 'KEscan'\
           else KExtract(directory, **kwargs)
