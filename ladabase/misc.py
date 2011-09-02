""" Miscellaeneous routinese for ladabase. """
def get_username():
  """ Returns username from $HOME/.lada file. """
  import lada

  if not hasattr(lada, "username"):
    raise RuntimeError("Cannot push OUTCAR if nobody is to blame.\n"\
                       "Please add 'username = \"your name\"' to $HOME/.lada.")
def get_ladabase(): 
  """ Return connector to the database. """
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()
  assert 'ladabase' in ip.user_ns, \
         RuntimeError("ladabase object not found in user namespace.")
  return ip.user_ns['ladabase']

