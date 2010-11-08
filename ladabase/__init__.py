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
    """ Pushes OUTCAR to database. 

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.
    """
    from os.path import exists
    from os.path import isfile
    from hashlib import sha512
    from ..opt import RelativeDirectory
    path = RelativeDirectory(filepath).path

    assert exists(path), ValueError('{0} does not exist.'.format(path))
    assert isfile(path), ValueError('{0} is not a file.'.format(path))

    with open(path, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()

    if self.outcars.exists(sha512=hash): 
      print "{0} already in database. Please use 'ladabase.update'.".format(filepath)
      return
    if 'filename' not in kwargs: kwargs['filename'] = filepath
    return self.outcars.put(string, sha512=hash, **kwargs)

  def update(self, filepath, **kwargs):
    """ Update a database entry. 

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.
    """
    from os.path import exists, isfile
    from hashlib import sha512
    from ..opt import RelativeDirectory

    path = RelativeDirectory(filepath).path
    assert exists(path), ValueError('{0} does not exist.'.format(path))
    assert isfile(path), ValueError('{0} is not a file.'.format(path))

    with open(path, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()

    if not self.outcars.exists(sha512=hash): 
      print "{0} not found in database. Please use 'ladabase.push'.".format(filepath)
      return 
    self.files.update({'sha512': hash}, {'$set': kwargs})
    return self.files.find_one(sha512=hash)['_id']

  def delete(self, filepath):
    """ Deletes an entry corresponding to the input file. """
    from os.path import exists, isfile
    from hashlib import sha512
    from ..opt import RelativeDirectory

    path = RelativeDirectory(filepath).path
    assert exists(path), ValueError('{0} does not exist.'.format(path))
    assert isfile(path), ValueError('{0} is not a file.'.format(path))

    with open(path, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()

    for u in self.database.files.find(sha512=hash): self.outcars.delete(u['_id'])

  def find_fromfile(self, path):
    """ Returns the database object corresponding to this file.

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.

        Finds the corresponding file using sha512 hash. 
    """
    from os.path import exists, isfile
    from hashlib import sha512
    from ..opt import RelativeDirectory

    ipath = RelativeDirectory(path).path
    assert exists(ipath), IOError('{0} does not exist.'.format(path))
    assert isfile(ipath), IOError('{0} is not a file.'.format(path))

    with open(ipath, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()
   
    assert self.outcars.exists(sha512=hash),\
           ValueError('{0} could not be found in database.'.format(path))

    return self.files.find_one(sha512=hash)


# class MassExtract(AbstractMassExtract):
#   """ Extracts a database of OUTCAR files. """
#   def __init__(self, path=None, comm=None, **kwargs):

def ipy_init():
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()

  manager = Manager()
  ip.user_ns['ladabase'] = manager
