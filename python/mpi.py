""" Package to handle mpi paradigma


    At present, we have mpi and non-mpi paradigma.
"""
# __docformat__ = "restructuredtext en"
__all__ = [ 'Communicator', 'world', 'all_gather', 'all_reduce', 'all_to_all',
            'broadcast', 'reduce', 'gather', 'NullComm']
from . import lada_with_mpi

def null_broadcast(comm, value=None, root=0): return value
def null_gather(comm, value, root=0): return [value]
def null_reduce(comm, value, op, root=0): return value
def null_split(comm, arg): return comm
def null_all_gather(comm, value): return [value]
def null_all_reduce(comm, value, op): return value
def null_all_to_all(comm, value=None): return [value]
def null_barrier(comm): pass


if not lada_with_mpi:
  broadcast = null_broadcast
  gather = null_gather
  reduce = null_reduce
  split = null_split
  all_gather = null_all_gather
  all_reduce = null_all_reduce
  all_to_all = null_all_to_all
else: 
  from boost.mpi import Communicator as BoostComm, world, all_gather, all_reduce, all_to_all,\
                        broadcast as boost_broadcast, reduce as boost_reduce,\
                        gather as boost_gather

  def broadcast(comm, value=None, root=0):
    return boost_broadcast(comm, value, root) if comm.is_mpi else null_broadcast(comm, value)
  broadcast.__doc__ = boost_broadcast.__doc__
  def reduce(comm, value, op, root=0):
    return boost_reduce(comm, value, op, root) if comm.is_mpi else null_reduce(comm, value, op)
  reduce.__doc__ = boost_reduce.__doc__
  def gather(comm, value, root=0):
    return boost_gather(comm, value, root) if comm.is_mpi else null_gather(comm, value)
  gather.__doc__ = boost_gather.__doc__
  def real(self):
    """ True if a real mpi-communicator. """
    return True
  def is_root(self):
    """ True if self.rank == 0. """
    return self.rank == 0
  def is_mpi(self):
    """ True if more than one proc. """
    return self.size > 1
  BoostComm.all_gather = all_gather
  BoostComm.all_reduce = all_reduce
  BoostComm.all_to_all = all_to_all
  BoostComm.broadcast = broadcast
  BoostComm.reduce = reduce
  BoostComm.gather = gather
  BoostComm.real = property(real)
  BoostComm.is_root = property(is_root)
  BoostComm.is_mpi = property(is_mpi)

class NullComm(object):
  """ Fake communicator for non-mpi stuff. """
  @property
  def is_mpi(self): return False
  @property
  def is_root(self): return True
  @property
  def rank(self): return 0
  @property
  def size(self): return 1
  @property 
  def real(self): return False
  
  def __init__(self): object.__init__(self)
  broadcast  = null_broadcast
  gather     = null_gather
  reduce     = null_reduce
  split      = null_split
  all_gather = null_all_gather
  all_reduce = null_all_reduce
  all_to_all = null_all_to_all
  def barrier(self): pass
  def __eq__(self, value): return value == None or isinstance(value, NullComm)
  def __ne__(self, value): return not self.__eq__(value)
  def __setstate__(self, value): self.__dict__.update(value)
  def __getstate__(self): return self.__dict__.copy()


if not lada_with_mpi: world, Communicator = NullComm(), NullComm

def Communicator(comm = None, with_world=False): 
  """ Communication object.


      In most case, this will wrap either none or a boost.mpi communicator.
  """
  if comm == None: return world if with_world else NullComm() 
  if isinstance(comm, BoostComm): return comm
  if isinstance(comm, NullComm): return comm
  raise ValueError("Unknown communicator type.")
