""" PyMongo interface. """
__pymongo_host__ = 'localhost'
__pymongo_port__ = 27017
__vasp_database_name__ = 'vasp'
__OUTCARS_prefix__ = 'OUTCARs'

from ..jobs import AbstractMassExtract
from ..opt import AbstractExtractBase, OutcarSearchMixin
from ..vasp.extract._common import Extract as ExtractCommonBase
from ..vasp.extract._dft import Extract as ExtractDFTBase
from ..vasp.extract._gw import Extract as ExtractGWBase

def get_username():
  """ Returns username from $HOME/.lada file. """
  import lada

  if not hasattr(lada, "username"):
    raise RuntimeError("Cannot push OUTCAR if nobody is to blame.\n"\
                       "Please add 'username = \"your name\"' to $HOME/.lada.")
def get_ladabase(): 
  """ Return connector to the database. """
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()
  assert 'ladabase' in ip.user_ns, \
         RuntimeError("ladabase object not found in user namespace.")
  return ip.user_ns['ladabase']

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
  def funccars_prefix(self):
    """ Name of the OUTCAR GridFS collection. """
    return self._outcars_prefix

  @property 
  def files(self):
    """ Return the OUTCAR files collection. """
    return getattr(self.database, '{0}.files'.format(self.outcars_prefix))

  def push(self, filename, outcar, **kwargs):
    """ Pushes OUTCAR to database. 

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.
    """
    from hashlib import sha512
    from os import uname


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

class StringStream(object):
  """ Converts string to file-like object. """
  def __init__(self, *args, **kwargs):
    """ Initializes the object using the stream. """
    from cStringIO import StringIO
    super(StringStream, self).__init__()
    self._stringio = StringIO(*args, **kwargs)
    """ String-as-a-filestream. """

  def __enter__(self):
    """ Makes this a context returning a cStringIO. """
    return self._stringio
  def __exit__(self, *args, **kwargs):
    """ Makes this a context returning a cStringIO. """
    self._stringio.close()

  def __getattr__(self, name):
    """ Forwards to stringio. """
    return getattr(self._stringio, name)

class ExtractCommon(AbstractExtractBase, ExtractCommonBase, OutcarSearchMixin):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, filter, comm=None, **kwargs):
    """ Initializes extraction object. """
    AbstractExtractBase.__init__(self, comm=comm)
    ExtractCommonBase.__init__(self)
    OutcarSearchMixin.__init__(self)
    filter = filter if isinstance(filter, dict) else {'_id': filter}
    n = [u for u in get_ladabase().files.find(filter)]
    if len(n) != 1: raise ValueError("{0} OUTCARS found from current filter.\n".format(len(n)))
    self._element = n[0]
    """ Elements in database associated with current OUTCAR. """

  def __outcar__(self):
    if self.compression == "bz2":  
      from bz2 import decompress
      return StringStream(decompress(get_ladabase().outcars.get(self._element['_id']).read()))
    elif self.compression.lower() == "none":
      return get_ladabase().outcars.get(self._element['_id'])
  def __funccar__(self):
    raise IOError('FUNCCARs are not stored in the database.')
  def __contcar__(self):
    raise IOError('CONTCARs are not stored in the database.')

  def __getattr__(self, name):
    """ Adds read-only properties corresponding to database metadata. """
    if name in self._element: return self._element[name]
    raise AttributeError()

  def __dir__(self):
    """ Attributes in this object. """
    result = set([u for u in self._element.keys() if u[0] != '_'] + ['_id', 'success']) \
             | set([u for u in dir(self.__class__) if u[0] != '_'])
    return list(result)


  @property
  def success(self):
    """ True if calculation was successfull. """
    return ExtractCommonBase.success.__get__(self)

class ExtractDFT(ExtractCommon, ExtractDFTBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, filter, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, filter, comm, **kwargs)
    ExtractDFTBase.__init__(self)

class ExtractGW(ExtractCommon, ExtractGWBase):
  """ Extracts GW data from an OUTCAR. """
  def __init__(self, filter, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, filter, comm, **kwargs)
    ExtractGWBase.__init__(self)

def VaspExtract(filter, *args, **kwargs): 
  """ Chooses between DFT or GW extraction object, depending on OUTCAR. """
  a = ExtractCommon(filter, *args, **kwargs)
  try: which = ExtractDFT if a.is_dft else ExtractGW
  except: which = ExtractCommon
  return which(filter, *args, **kwargs)

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
    # this will throw on unknown kwargs arguments.
    super(MassExtract, self).__init__(**kwargs)

    self.Extract = Extract if Extract != None else VaspExtract
    """ Extraction class to use. """


  @property 
  def ladabase(self): return get_ladabase()

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
  from .ipython import push, amend
  ip = get_ipy()

  manager = Manager()
  ip.user_ns['ladabase'] = manager
  ip.expose_magic("push", push)
  ip.expose_magic("amend", amend)

