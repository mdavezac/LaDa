""" PyMongo interface. """
__pymongo_host__ = 'localhost'
__pymongo_port__ = 27017
__vasp_database_name__ = 'vasp'
__OUTCARS_collection_name__ = 'OUTCARs'

class Manager(object): 
  """ Holds data regarding database management. """
  def __init__(self, host=None, port=None, database=None, collection=None): 
    """ Initializes a connection and database. """
    from pymongo import Connection
    from gridfs import GridFs
    super(Manager, self).__init__()

    self._host = host if host != None else __pymongo_host__
    """ Host where the database is hosted. """
    self._port = port if port != None else __pymongo_port__
    """ Port of the host where the database is hosted. """
    self._vaspbase_name = database if database != None else __vasp_database_name__
    """ Name of the vasp database. """
    self._outcarbase_name = name if name != None else __OUTCARS_collection_name__
    """ Name of the OUTCAR database. """
    self.connection = Connection(self._host, self._port)
    """ Holds connection to pymongo. """
    self.database = getattr(self.connection, self._vaspbase_name)
    """ Database within pymongo. """
    self.outcars = GridFs(self.database, self._outcarbase_name)
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
  def outcar_collection_name(self):
    """ Name of the OUTCAR GridFS collection. """
    return self._outcarbase_name

def outcar_push(self, arg):
  """ Push outcar to database. """
  from hashlib import sha224
  
  if not 'ladabase_manager' in ip.user_ns:
    print "Could not find ladabase in user namespace."
    return
  ladabase = ip.user_ns['ladabase_manager']

  filename = args.rstrip().lstrip()
  with open(filename, 'r') as file: string = file.read()
  hash = sha224(string).hexdigest()

  if ladabase.outcars.exists(hash=hash): return
  return ladabase.outcars.put(string, filename=filename, hash=hash)


def ipy_init():
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()

  manager = Manager()
  ip.user_ns['ladabase_manager'] = manager
