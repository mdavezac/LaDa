""" Holds a list of callbacks to call before exit. 

    Helps declares functions which are called if python exists abruptly.
    These functions are given a unique identifier so they can be removed if not
    needed anymore.
"""
from atexit import register

_list_of_callbacks = {}
""" List of callbacks + arguments. """


@register
def _atexit_onexit():
  """ Specific at-exit function for lada. """
  global _list_of_callbacks
  _list_of_callbacks.pop('abort', None)
  _list_of_callbacks.pop('term', None)
  for callback, args, kwargs in _list_of_callbacks.values():
    if callback is not None: 
      try: callback(*args, **kwargs)
      except: pass

def _onexit_signal(signum, stackframe):
  from signal import SIGABRT, SIGTERM, signal, SIG_DFL
  abort = _list_of_callbacks.get('abort', None)
  term = _list_of_callbacks.get('term', None)
  if signum == SIGABRT and abort is not None:
    try: signal(SIGABRT, abort)
    except: signal(SIGABRT, SIG_DFL)
  elif signum == SIGTERM and term is not None: 
    try: signal(SIGTERM, term)
    except: signal(SIGTERM, SIG_DFL)
  else: signal(SIGABRT, SIG_DFL)
  raise SystemExit(signum)

# delete register from atexit
del register

def add_callback(callback, *args, **kwargs):
  """ Adds function to call at exit. 
  
      :param callback:
        Callback function.
      :param *args: 
        Arguments to the callback.
      :param **kwargs:
        Keyword arguments to the callback.

      :return: An identifier with which the callback can be deleted.
  """
  from uuid import uuid4
  id = uuid4()
  _list_of_callbacks[id] = callback, args, kwargs
  return id

def del_callback(id):
  """ Deletes a callback from the list. """
  _list_of_callbacks.pop(id, None)

# on first opening this module, change sigterm signal.
if len(_list_of_callbacks) == 0:
  from signal import SIGABRT, SIGTERM, signal
  _list_of_callbacks['abort'] = signal(SIGABRT, _onexit_signal)
  _list_of_callbacks['term'] = signal(SIGTERM, _onexit_signal)
  del signal
  del SIGABRT
  del SIGTERM
