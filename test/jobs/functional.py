# Functional which will be called
class Functional(object):
  def __call__(self, *args, **kwargs):
    from time import sleep
    from boost.mpi import world
    
    if "wait" in kwargs: sleep(kwargs["wait"])
    comm = kwargs.pop("comm", None)
    return "Rank: %i\nargs: %s\nkwargs: %s" % (world.rank, args, kwargs)


