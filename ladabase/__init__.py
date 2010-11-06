""" PyMongo interface. """
__pymongo_host__ = 'localhost'
__pymongo_port__ = 27017
__vasp_database_name__ = 'vasp'
__OUTCARS_prefix__ = 'OUTCARs'

class Manager(object): 
  """ Holds data regarding database management. """
  def __init__(self, host=None, port=None, database=None, prefix=None): 
    """ Initializes a connection and database. """
    from pymongo import Connection
    from gridfs import GridFS
    super(Manager, self).__init__()

    self._host = host if host != None else __pymongo_host__
    """ Host where the database is hosted. """
    self._port = port if port != None else __pymongo_port__
    """ Port of the host where the database is hosted. """
    self._vaspbase_name = database if database != None else __vasp_database_name__
    """ Name of the vasp database. """
    self._outcars_prefix = prefix if prefix != None else __OUTCARS_prefix__
    """ Name of the OUTCAR database. """
    self.connection = Connection(self._host, self._port)
    """ Holds connection to pymongo. """
    self.database = getattr(self.connection, self._vaspbase_name)
    """ Database within pymongo. """
    self.outcars = GridFS(self.database, self.outcars_prefix)
    """ GridFS object for OUTCARs. """

  @property
  def host(self):
    """ Host where the database is hosted. """
    return self._host
  @property
  def port(self):
    """ Port of the host where the database is hosted. """
    return self._port
  @property
  def vasp_database_name(self):
    """ Name of the vasp database. """
    return self._vaspbase_name
  @property
  def outcars_prefix(self):
    """ Name of the OUTCAR GridFS collection. """
    return self._outcars_prefix

  @property 
  def files(self):
    """ Return the files collection. """
    return getattr(self.database, '{0}.files'.format(self.outcars_prefix))

  def push(self, filepath, **kwargs):
    """ Pushes OUTCAR to database. """
    from os.path import exists
    from hashlib import md5
    from ..opt import RelativeDirectory
    path = RelativeDirectory(filepath).path

    assert exists(path), ValueError('{0} does not exist.'.format(path))

    with open(path, 'r') as file: string = file.read()
    hash = md5(string).hexdigest()

    if 'filename' not in kwargs: kwargs['filename'] = filepath
    if self.outcars.exists(md5=hash): 
      self.files.update({'md': hash}, {'$set': kwargs})
      return self.files.find_one(md5=hash)['_id']
    return self.outcars.put(string, **kwargs)

  def find_fromfile(self, path):
    """ Returns the database object corresponding to this file.

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.

        Finds the corresponding file using md5 hash. 
    """
    from hashlib import md5
    from os.path import exists, isfile
    from ..opt import RelativeDirectory

    ipath = RelativeDirectory(path).path
    assert exists(ipath), IOError('{0} does not exist.'.format(path))
    assert isfile(ipath), IOError('{0} is not a file.'.format(path))

    with open(ipath, 'r') as file: string = file.read()
    hash = md5(string).hexdigest()
   
    assert self.outcars.exists(md5=hash),\
           ValueError('{0} could not be found in database.'.format(path))

    return self.files.find_one(md5=hash)

def ipy_init():
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()

  manager = Manager()
  ip.user_ns['ladabase'] = manager
