""" PyMongo interface. """
__pymongo_host__ = 'localhost'
__pymongo_port__ = 27017
__vasp_database_name__ = 'vasp'
__OUTCARS_prefix__ = 'OUTCARs'

from ..jobs import AbstractMassExtract
from ..opt import AbstractExtractBase
from ..vasp.extract._common import Extract as ExtractCommonBase
from ..vasp.extract._dft import Extract as ExtractDFTBase
from ..vasp.extract._gw import Extract as ExtractGWBase
from ..vasp.extract._mixin import SearchMixin


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
    from getpass import getuser
    from os import uname
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
    if 'uploader' not in kwargs: kwargs['uploader'] = getuser()
    if 'host'     not in kwargs: kwargs['host']     = uname()[1]
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
    return self.files.find_one({'sha512':hash})['_id']

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

    for u in self.database.files.find({'sha512':hash}): self.outcars.delete(u['_id'])

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

class StreamFSDoc(object):
  """ Read-only streamed database file which behaves like  a file object. """
  def __init__(self, _id):
    """ Initializes the object using the stream. """
    from cStringIO import StringIO
    object.__init__(self)
    self._id = _id['_id'] if isinstance(_id, dict) else _id
    a = [f for f in self.ladabase.files.find({'_id': self._id})] 
    assert len(a) == 1, RuntimeError('Could not find object')
    self._stringio = StringIO(self.ladabase.outcars.get(self._id).read())

  def __enter__(self):
    """ Makes this a context returning a cStringIO. """
    return self._stringio
  def __exit__(self, *args, **kwargs):
    """ Makes this a context returning a cStringIO. """
    self._stringio.close()

  def __getattr__(self, name):
    """ Forwards to stringio. """
    return getattr(self.stringio, name)

  @property 
  def ladabase(self):
    from IPython.ipapi import get as get_ipy
    ip = get_ipy()
    assert 'ladabase' in ip.user_ns, \
           RuntimeError("ladabase object not found in user namespace.")
    return ip.user_ns['ladabase']

class IOMixin(object): 
  """ Defines IO for ladabase object. """
  def __init__(self, id_or_doc):
    object.__init__(self)
    self._id = id_or_doc['_id'] if isinstance(id_or_doc, dict) else id_or_dcoc
    
  def __outcar__(self):
    return StreamFSDoc(self._id)
  def __funccar__(self):
    raise IOError('FUNCCARs are not stored in the database.')
  def __contcar__(self):
    raise IOError('CONTCARs are not stored in the database.')


class ExtractCommon(AbstractExtractBase, ExtractCommonBase, IOMixin, SearchMixin):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, id_or_doc, comm=None, **kwargs):
    """ Initializes extraction object. """
    from os.path import exists, isdir, basename, dirname
    AbstractExtractBase.__init__(self, comm=comm)
    ExtractCommonBase.__init__(self)
    IOMixin.__init__(self, id_or_doc, **kwargs)
    SearchMixin.__init__(self)

class ExtractDFT(ExtractCommon, ExtractDFTBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, id_or_doc, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, id_or_doc, comm, **kwargs)
    ExtractDFTBase.__init__(self)

class ExtractGW(ExtractCommon, ExtractGWBase):
  """ Extracts GW data from an OUTCAR. """
  def __init__(self, id_or_doc, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, id_or_doc, comm, **kwargs)
    ExtractGWBase.__init__(self)

def VaspExtract(id_or_doc, *args, **kwargs): 
  """ Chooses between DFT or GW extraction object, depending on OUTCAR. """
  a = ExtractCommon(id_or_doc, *args, **kwargs)
  try: which = ExtractDFT if a.is_dft else ExtractGW
  except: which = ExtractCommon
  return which(id_or_doc, *args, **kwargs)

class MassExtract(AbstractMassExtract):
  """ Propagates vasp extractors from all subdirectories.
  
      Trolls through all subdirectories for vasp calculations, and organises
      results as a dictionary where keys are the name of the diretory.
  """
  def __init__(self, Extract = None, **kwargs):
    """ Initializes MassExtract.
    
    
        :Parameters:
          Extract : `lada.vasp.Extract`
            Extraction class to use. 
          kwargs : dict
            Keyword parameters passed on to AbstractMassExtract.

        :kwarg naked_end: True if should return value rather than dict when only one item.
        :kwarg unix_re: converts regex patterns from unix-like expression.
    """
    from os import getcwd
    from os.path import exists, isdir
    from ..opt import RelativeDirectory

    # this will throw on unknown kwargs arguments.
    super(MassExtract, self).__init__(**kwargs)

    self.Extract = Extract if Extract != None else VaspExtract
    """ Extraction class to use. """


  @property 
  def ladabase(self):
    from IPython.ipapi import get as get_ipy
    ip = get_ipy()
    assert 'ladabase' in ip.user_ns, \
           RuntimeError("ladabase object not found in user namespace.")
    return ip.user_ns['ladabase']

  def __iter_alljobs__(self):
    """ Goes through all directories with an OUTVAR. """

    for doc in self.ladabase.files.find():
      yield str(doc['_id']), self.Extract(doc)

  def __copy__(self):
    """ Returns a shallow copy. """
    result = self.__class__(self.Extract)
    result.__dict__.update(self.__dict__)
    return result


def ipy_init():
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()

  manager = Manager()
  ip.user_ns['ladabase'] = manager
