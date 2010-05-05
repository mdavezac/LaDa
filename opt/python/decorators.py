""" A subpackage with decorators """


def broadcast_result(key=False, attr=False, which=0): 
  """ Decorator to wrap method such that its result is broadcasted to all nodes. """
  def _key_broadcast_result(method):
    """ Decorator to wrap method such that its result is broadcasted to all nodes.
         
        The method is performed on on the root node and the results broadcasted to all.
        The returned method accepts a communicator via the keyword argument comm.
        If comm is not present or comm is None, then all processes call the
        original method, and the results are not broadcasted.
    """
    from boost.mpi import world, broadcast

    def wrapped(*args, **kwargs):
      # no comm provided, performs on all procs
      if "comm" not in kwargs: return method(*args, **kwargs)
      # removes comm argument
      comm = kwargs["comm"]
      del kwargs["comm"]
      # nonetype case: each proc performs same action. 
      if comm == None: return method(*args, **kwargs)
      # not an mpi process.
      if comm.size == 1: return method(*args, **kwargs)
      # is an mpi process.
      error, result, exception = False, None, None
      if comm.rank == 0: # root process
        try: result = method(*args, **kwargs)
        except Exception as exception: error, result = True, (str(exception), world.rank)
        broadcast(comm, (error, result), root=0)
        assert not error, exception
      else:
        error, result = broadcast(comm, root=0)
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
    from boost.mpi import world, broadcast

    def wrapped(*args, **kwargs):
      assert len(args) > which,\
             RuntimeError("Expected at least %i arguments, got %s." % (which, args))
      assert hasattr(args[which], "comm"),\
             RuntimeError("Argument %i does not have communicator." %(which))
      # nonetype case: each proc performs same action. 
      if args[which].comm == None: return method(*args, **kwargs)
      # not an mpi process.
      if args[which].comm.size == 1: return method(*args, **kwargs)
      # is an mpi process.
      error, result, exception = False, None, None
      if args[which].comm.rank == 0: # root process
        try: result = method(*args, **kwargs)
        except Exception as exception: error, result = True, (str(exception), world.rank)
        broadcast(args[which].comm, (error, result), root=0)
        assert not error, exception
      else:
        error, result = broadcast(args[which].comm, root=0)
        assert not error, RuntimeError("Process %i reports an error: %s"  % (result[1], result[0]))
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

def add_setter(method, docstring): 
  """ Adds an input-like setter property. """
  def _not_available(self): raise RuntimeError("Error: No cheese available.")
  return property(_not_available, method, docstring)
