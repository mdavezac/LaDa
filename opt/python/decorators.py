""" A subpackage with decorators """
import boost.mpi

def broadcast_result(method):
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
  def wrapped(*args, **kwargs):
    # no comm provided, performs on all procs
    if "comm" not in kwargs: return method(*args, **kwargs)
    # removes comm argument
    comm = kwargs["comm"]
    del kwargs["comm"]
    # not an mpi process.
    if comm.size == 1: return method(*args, **kwargs)
    # is an mpi process.
    error, result = False, None
    if comm.rank == 0: # root process
      try: result = method(*args, **kwargs)
      except Exception as inst: error, result = True, (inst, boost.mpi.world.rank)
      boost.mpi.broadcast(comm, (error, result), root=0)
      if error: raise result
    else:
      error, result = boost.mpi.broadcast(comm, root=0)
      if error:
        raise RuntimeError,\
              "Process %i reports an error: %s"  % (result[1], result[0])
    return result
  wrapped.__name__ = method.__name__
  wrapped.__doc__ = method.__doc__
  wrapped.__module__ = method.__module__
  wrapped.__dict__.update(method.__dict__)
  return wrapped
 
def make_cached(method):
  """ Caches the result of a method. """

  def wrapped(*args, **kwargs):
    assert len(args) > 0, "expected bound method."
    cache_name = "_cached_attr%s" % (method.__name__)
    if not hasattr(args[0], cache_name): 
      setattr(args[0], cache_name, method(*args, **kwargs))
    return getattr(args[0], cache_name)

  wrapped.__name__ = method.__name__
  wrapped.__doc__ = method.__doc__
  wrapped.__module__ = method.__module__
  wrapped.__dict__.update(method.__dict__)
  return wrapped

def uncache(ob):
  """ Uncaches results cached by @make_cached. """ 
  for key in ob.__dict__.keys():
    if key[:len("_cached_attr")] == "_cached_attr": del ob.__dict__[key]
