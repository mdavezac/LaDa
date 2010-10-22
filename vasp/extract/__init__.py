""" Subpackage containing extraction methods for vasp parameters from vasp output. """
__docformat__  = 'restructuredtext en'
from ._dft import Extract
from ._gw import  ExtractGW
__all__ = ['Extract', 'ExtractGW']

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

          :Kwargs naked_end: True if should return value rather than dict when only one item.
          :Kwargs unix_re: converts regex patterns from unix-like expression.
      """
      from os import getcwd
      from os.path import exists, isdir
      from . import Extract as VaspExtract
      from ...opt import RelativeDirectory

      # this will throw on unknown kwargs arguments.
      super(MassExtract, self).__init__(**kwargs)

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

    def walk_through(self):
      """ Goes through all directories with a contcar. """
      from os import walk, getcwd
      from os.path import abspath, relpath, abspath, join

      for dirpath, dirnames, filenames in walk(self.rootdir, topdown=True, followlinks=True):
        if self.OUTCAR not in filenames: continue

        try: result = self.Extract(join(self.rootdir, dirpath))
        except: continue

        result.OUTCAR = self.OUTCAR
        yield join('/', relpath(dirpath, self.rootdir)), result

    @property
    def _attributes(self): 
      """ Returns __dir__ set special to the extraction itself. """
      return set([u for u in dir(self.Extract()) if u[0] != '_'])

    def __copy__(self):
      """ Returns a shallow copy. """
      result = self.__class__(self.rootdir)
      result.__dict__.update(self.__dict__)
      return result

      


