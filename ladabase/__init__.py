""" PyMongo interface. """
__all__ = ['Manager', 'MassExtract', 'VaspExtract', 'ipy_init']

from .vasp import VaspExtract
from .massextract import MassExtract

class Manager(object): 
  """ Holds data regarding database management. """
  def __init__(self, host=None, port=None, database=None, prefix=None): 
    """ Initializes a connection and database. """
    from pymongo import Connection
    from gridfs import GridFS
    from .. import pymongo_host, pymongo_port, vasp_database_name,  OUTCARS_prefix
    super(Manager, self).__init__()

    self._host = host if host != None else pymongo_host
    """ Host where the database is hosted. """
    self._port = port if port != None else pymongo_port
    """ Port of the host where the database is hosted. """
    self._vaspbase_name = database if database != None else vasp_database_name
    """ Name of the vasp database. """
    self._outcars_prefix = prefix if prefix != None else OUTCARS_prefix
    """ Name of the OUTCAR database. """
    self.connection = Connection(self._host, self._port)
    """ Holds connection to pymongo. """
    self.database = getattr(self.connection, self._vaspbase_name)
    """ Database within pymongo. """
    self.outcars = GridFS(self.database, self.outcars_prefix)
    """ GridFS object for OUTCARs. """
    self.comments = self.database["{0}.comments".format(self.outcars_prefix)]
    """ Collection of comments attached when adding OUTCAR's to a file. """
    self.files = self.database["{0}.files".format(self.outcars_prefix)]
    """ OUTCAR files collection. """
    self.extracted = self.database["extracted"]
    """ Collection with pre-extracted values from the outcar. """

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

  def push(self, filename, outcar, comment, **kwargs):
    """ Pushes OUTCAR to database. 

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.
    """
    from hashlib import sha512
    from os import uname
    from .misc import get_username

    assert len(comment.replace(' ', '').replace('\n', '')) != 0,\
           ValueError("Cannot push file if comment is empty.")
    
    try: kwargs["comment"] = self.comments.find({'text': comment}).next()
    except StopIteration: # add comment to database
      kwargs["comment"] = self.comments.insert({'text': comment})
      print kwargs['comment']

    hash = sha512(outcar).hexdigest()
    if self.outcars.exists(sha512=hash): 
      print "{0} already in database. Please use 'ladabase.update'.".format(filename)
      return 
    if 'filename' not in kwargs: kwargs['filename'] = filename
    if 'uploader' not in kwargs: kwargs['uploader'] = get_username()
    if 'host'     not in kwargs: kwargs['host']     = uname()[1]
    compression = kwargs.get('compression', None)
    kwargs['compression'] = compression
    if compression == "bz2": 
      from bz2 import compress
      return self.outcars.put(compress(outcar), sha512=hash, **kwargs)
    elif compression == None: return self.outcars.put(outcar, sha512=hash, **kwargs)
    else: raise ValueError("Invalid compression format {0}.".format(compression))

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

    return self.files.find_one({'sha512': hash})

  def __contains__(self, path):
    """ True if path already in database. """
    from os.path import exists
    from hashlib import sha512
    from ..opt import RelativeDirectory
    path = RelativeDirectory(path).path
    if not exists(path): ValueError("File {0} does not exist.".format(path))
    with open(path, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()
    return self.outcars.exists(sha512=hash)



def ipy_init():
  from IPython.ipapi import get as get_ipy
  from .ipython import push, amend
  from .filter import init as filter_init
  ip = get_ipy()

  manager = Manager()
  ip.user_ns['ladabase'] = manager
  ip.expose_magic("push", push)
  ip.expose_magic("amend", amend)
  filter_init(ip)

