""" Holds a list of callbacks to call before exit. 

    Helps declares functions which are called if python exists abruptly.
    These functions are given a unique identifier so they can be removed if not
    needed anymore.
"""
from atexit import register

_callback_dict = {}
""" List of callbacks + arguments. """


def _call_callbacks():
  """ Calls on-exit functions. """
  try: 
    _callback_dict.pop('abort', None)
    _callback_dict.pop('term', None)
    while len(_callback_dict):
      name, (callback, args, kwargs) = _callback_dict.popitem()
      if callback is not None: 
        try: callback(*args, **kwargs)
        except: pass
  except: pass

@register
def _atexit_onexit(): 
  """ Specific at-exit function for pylada. """
  _call_callbacks()


def _onexit_signal(signum, stackframe):
  from signal import SIGABRT, SIGTERM, signal, SIG_DFL

  abort = _callback_dict.pop('abort', None)
  term  = _callback_dict.pop('term', None)

  _call_callbacks()

  if signum == SIGABRT and abort is not None:
    try: signal(SIGABRT, abort)
    except: signal(SIGABRT, SIG_DFL)
  elif signum == SIGTERM and term is not None: 
    try: signal(SIGTERM, term)
    except: signal(SIGTERM, SIG_DFL)
  else: signal(SIGABRT, SIG_DFL)
  raise SystemExit(signum)

# delete register from module
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
  _callback_dict[id] = callback, args, kwargs
  return id

def del_callback(id):
  """ Deletes a callback from the list. """
  _callback_dict.pop(id, None)

# on first opening this module, change sigterm signal.
if len(_callback_dict) == 0:
  from signal import SIGABRT, SIGTERM, signal
  _callback_dict['abort'] = signal(SIGABRT, _onexit_signal)
  _callback_dict['term'] = signal(SIGTERM, _onexit_signal)
  del signal
  del SIGABRT
  del SIGTERM
