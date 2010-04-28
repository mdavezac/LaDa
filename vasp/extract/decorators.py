""" General decorators to use in extraction routines. """

def bound_broadcast_result(method):
  """ Decorator to wrap method such that its result is broadcasted to all nodes.
  
      The method should be a bound method. Some constraints are imposed on the
      method. Use with care.
      The method is performed on on the root node and the results broadcasted to all.
      The returned method accepts a communicator via the keyword argument comm.
      If comm is not present or comm is None, then all processes call the
      original method, and the results are not broadcasted.
  """
  def wrapped(*args, **kwargs):
    from boost.mpi import broadcast, world
    assert len(args) > 0, RuntimeError("Expected a bound method.")
    # no comm provided, performs on all procs
    if not hasattr(args[0],"comm"): return method(*args, **kwargs)
    # not an mpi process.
    comm = args[0].comm
    if comm == None: return method(*args, **kwargs)
    if comm.size == 1: return method(*args, **kwargs)
    # is an mpi process.
    error, result = False, None
    if comm.rank == 0: # root process
      try: result = method(*args, **kwargs)
      except Exception as inst: error, result = True, (inst, world.rank)
      broadcast(comm, (error, result), root=0)
      if error: raise result
    else:
      error, result = broadcast(comm, root=0)
      if error:
        raise RuntimeError,\
              "Process %i reports an error: %s"  % (result[1], result[0])
    return result
  wrapped.__name__ = method.__name__
  wrapped.__doc__ = method.__doc__
  wrapped.__module__ = method.__module__
  wrapped.__dict__.update(method.__dict__)
  return wrapped
