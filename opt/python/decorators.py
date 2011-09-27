""" A subpackage with decorators """
__docformat__ = "restructuredtext en"
__all__ = ['FileCache', 'broadcast_result', 'count_calls', 'make_cached', 'add_cache', 'add_setter']
from .filecache import FileCache


def broadcast_result(key=False, attr=False, which=0): 
  """ Decorator to wrap method such that its result is broadcasted to all nodes. """
  def _key_broadcast_result(method):
    """ Decorator to wrap method such that its result is broadcasted to all nodes.
         
        The method is performed on on the root node and the results broadcasted to all.
        The returned method accepts a communicator via the keyword argument comm.
        If comm is not present or comm is None, then all processes call the
        original method, and the results are not broadcasted.
    """

    def wrapped(*args, **kwargs):
      from ..mpi import Communicator, world
      # removes comm argument
      comm = Communicator(kwargs.pop("comm", None))
      # nonetype case: each proc performs same action. Serial case.
      if not comm.is_mpi: return method(*args, **kwargs)
      # is an mpi process.
      error, result, exception = False, None, None
      if comm.is_root:
        try: result = method(*args, **kwargs)
        except Exception as exception:
          comm.broadcast((True, (str(exception), world.rank)))
          raise
        comm.broadcast((False, result))
      else:
        error, result = comm.broadcast()
        assert not error, RuntimeError("Process %i reports an error: %s"  % (result[1], result[0]))
      return result
    wrapped.__name__ = method.__name__
    wrapped.__doc__ = method.__doc__
    wrapped.__module__ = method.__module__
    wrapped.__dict__.update(method.__dict__)
    return wrapped

  def _attr_broadcast_result(method):
    """ Decorator to wrap method such that its result is broadcasted to all nodes.
         
        The method is performed on on the root node and the results broadcasted to all.
        It is expected that a communicator name comm is found on the first
        argument of the method (eg in self).
    """

    def wrapped(*args, **kwargs):
      from ..mpi import Communicator, world
      assert len(args) > which,\
             RuntimeError("Expected at least %i arguments, got %s." % (which, args))
      assert hasattr(args[which], "comm"),\
             RuntimeError("Argument %i does not have communicator." %(which))
      comm = Communicator(args[which].comm)
      # nonetype and serial case: each proc performs same action. 
      if not comm.is_mpi: return method(*args, **kwargs)
      # is an mpi process.
      error, result, exception = False, None, None
      if comm.is_root:
        try: result = method(*args, **kwargs)
        except Exception as exception: error, result = True, (str(exception), world.rank)
        comm.broadcast((error, result))
        assert not error, exception
      else:
        error, result = comm.broadcast(comm)
        assert not error, RuntimeError("Process %i reports an error: %s"  % (result[1], result[0]))
      if __debug__: comm.barrier()
      return result
    wrapped.__name__ = method.__name__
    wrapped.__doc__ = method.__doc__
    wrapped.__module__ = method.__module__
    wrapped.__dict__.update(method.__dict__)
    return wrapped

  assert (key or attr) and (not (key and attr)),\
         RuntimeError("comm in broadcast_result should be bound either to key or attribute.")
  if key: return _key_broadcast_result
  return _attr_broadcast_result
 
def count_calls(name='nbcalc', argnb=0):
  """ Counts calls to a member function. """
  def decorating_function(method):

    def wrapped(*args, **kwargs):
      if not hasattr(args[argnb], name): setattr(args[argnb], name, 0)
      result = method(*args, **kwargs)
      setattr(args[argnb], name, getattr(args[argnb], name)+1)
      return result
   
    wrapped.__name__ = method.__name__
    wrapped.__doc__ = method.__doc__
    wrapped.__module__ = method.__module__
    wrapped.__dict__.update(method.__dict__)
    return wrapped
  return decorating_function
 
def make_cached(method):
  """ Caches the result of a method for futur calls. """

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

def add_setter(method, docstring = None): 
  """ Adds an input-like setter property. """
  def _not_available(self): raise RuntimeError("Error: No cheese available.")
  if docstring == None and hasattr(method, "__doc__"): docstring = method.__doc__
  return property(fget=_not_available, fset=method,  doc=docstring)
