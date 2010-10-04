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
    def __init__(self, path, hook = None, comm = None, Extract = None):
      """ Initializes MassExtract.
      
      
          :Parameters:
            - `path` : root directory for which to investigate all subdirectories.
            - `hook` : Optional callable taking a path as an argument and
              returning True or False. If it returns False, that particular
              directory and its subdirectories will not be listed as a vasp
              directory.
            - `comm` : boost.mpi.communicator. Can be None.
            - Extract : Extraction class to use. Defaults to `lada.vasp.Extract.`
      """
      from os.path import exists, isdir
      from . import Extract as VaspExtract
      self.hook = hook
      """ Optional callable to exclude directories from extraction. 
      
          Callable takes a path as an argument and returns True or False. If it
          returns False, that particular ill not be listed as a vasp directory.
      """
      self.Extract = Extract if Extract != None else VaspExtract
      """ Extraction class to use. """

      super(MassExtract, self).__init__(path, comm = None)
      assert exists(self.root), RuntimeError("Path {0} does not exist.".format(self.root))
      assert isdir(self.root), RuntimeError("Path {0} is not a directory.".format(self.root))

    def walk_through(self):
      """ Goes through all directories with a contcar. """
      from os import walk, getcwd
      from os.path import abspath, relpath, abspath, join

      for dirpath, dirnames, filenames in walk(self.root, topdown=True, followlinks=True):
        if "OUTCAR" not in filenames: continue
        if hasattr(self.hook, "__call__") and not self.hook(join(self.root, dirpath)):
          while len(dirnames): dirnames.pop()
          continue

        try: result = self.Extract(join(self.root, dirpath), comm = self.comm)
        except: pass
        else: yield relpath(dirpath, self.root), result

    def _properties(self): 
      """ Returns cached __dir__ result. """
      return set([u for u in dir(self.Extract()) if u[0] != '_'])

      


