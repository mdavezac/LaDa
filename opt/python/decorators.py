""" A subpackage with decorators """
import boost.mpi
from decorator import decorator as _dec_dec

@_dec_dec
def broadcast_result(method, *args, **kwargs):
  """ Performs operation on one proc and broadcasts to others. 
       
      method is called on the root processor only. Any exception raised there
      will be propagated to other processors.
      @param method: The method to decorate.
      @param args: The arguments to the method.
      @param **kwargs: keyword arguments to the method. A comm keyword will be
         extracted from the list if it exists. It is expected to be an mpi
         communicator and will not be passed on to the method. If no comm
         keyword is found, then boost.mpi.communicator is used instead.
  """
  # extracts comm keyword from dictionary
  comm = boost.mpi.world
  if "comm" in kwargs: 
    comm = kwargs["comm"]
    del kwargs["comm"]
  # not an mpi process.
  if comm.size == 1: return method(*args, **kwargs)
  # is an mpi process.
  error, result = False, None
  if comm.rank == 0: # root process
    try: result = method(*args, **kwargs)
    except Exception as inst: error, result = True, inst
    boost.mpi.broadcast(comm, (error, result), root=0)
  else: error, result = boost.mpi.broadcast(comm, root=0)
  if error: raise result
  return result
  
