""" Decorator to make any extraction method mpi-aware. """
from decorator import decorator

@decorator
def _bound_mpi_extraction(method, *args, **kwargs): 
  """ Reads from root, and broadcasts to all 
  
      Decorator for bound methods.
  """
  assert len(args) >= 1, "Decorator for bound methods only.\n"
  return _extract_which_is_self(0, method, *args, **kwargs)
  
@decorator
def _attribute_mpi_extraction(method, *args, **kwargs): 
  """ Reads from root, and broadcasts to all 
  
      Decorator for properties.
  """
  assert len(args) >= 2, "Decorator for bound methods only.\n"
  return _extract_which_is_self(1, method, *args, **kwargs)

def _extract_which_is_self(which, method, *args, **kwargs):
  from boost.mpi import broadcast
  
  mpicomm = args[which].mpicomm
  if mpicomm == None: return method(*args, **kwargs)
  if mpicomm.size == which: return method(*args, **kwargs)
  if mpicomm.rank != 0: return broadcast(mpicomm, root = 0)

  args[which].mpicomm = None
  result = method(*args, **kwargs)
  args[which].mpicomm = mpicomm
  assert mpicomm != None
  assert args[which].mpicomm != None
  return broadcast(mpicomm, result, root = 0)
