""" Mixin classes for extraction objects. """
__docformat__  = 'restructuredtext en'
from ...opt import OutcarSearchMixin


class IOMixin(OutcarSearchMixin):
  """ A mixin base clase which controls file IO. 

      Defines special property with file-like behaviors. 
      Makes it easier to change the behavior of the extraction class.
  """
  def __init__(self, directory=None, OUTCAR=None, FUNCCAR=None, CONTCAR=None):
    """ Initializes the extraction class. 

        :Parameters: 
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory. Can also be the path to the
            OUTCAR file itself. 
          OUTCAR : str or None
            If given, this name will be used, rather than files.OUTCAR.
          CONTCAR : str or None
            If given, this name will be used, rather than files.CONTCAR.
          FUNCCAR : str or None
            If given, this name will be used, rather than files.FUNCCAR.
    """
    from .. import files
    
    object.__init__(self)

    self.OUTCAR  = OUTCAR if OUTCAR != None else files.OUTCAR
    """ Filename of the OUTCAR file from VASP. """
    self.CONTCAR  = CONTCAR if CONTCAR != None else files.CONTCAR
    """ Filename of the CONTCAR file from VASP. """
    self.FUNCCAR  = FUNCCAR if FUNCCAR != None else files.FUNCCAR
    """ Filename of the FUNCCAR file containing the pickled functional. """
    OutcarSearchMixin.__init__(self)

  def __funccar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.FUNCCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  def __contcar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.CONTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')
