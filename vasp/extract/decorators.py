""" General decorators to use in extraction routines. """

def bound_broadcast_result(method):
  """ Performs operation on root node, then broadcast results. 
       
      The method will be called  on the root process only. Any exception raised
      there will be propagated to other processes.
      @param self: The bound method to make mpi-aware. The first argurment to
         the bound method(self) should have a "comm" attribute. If it doesn't,
         then the method is called for each process.
      @param args: The arguments to the method.
      @param **kwargs: keyword arguments to the method. A comm keyword will be
         extracted from the list if it exists. It is expected to be an mpi
         communicator and will not be passed on to the method. If no comm
         keyword is found, then boost.mpi.communicator is used instead.
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
