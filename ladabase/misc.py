""" Miscellaeneous routinese for ladabase. """
def get_username():
  """ Returns username from $HOME/.pylada file. """
  import pylada

  if not hasattr(pylada, "username"):
    raise RuntimeError("Cannot push OUTCAR if nobody is to blame.\n"\
                       "Please add 'username = \"your name\"' to $HOME/.pylada.")
  return pylada.username

def get_ladabase(): 
  """ Return connector to the database. """
  from IPython.ipapi import get as get_ipy
  ip = get_ipy()
  if 'ladabase' not in ip.user_ns: 
     raise RuntimeError("ladabase object not found in user namespace.")
  return ip.user_ns['ladabase']

