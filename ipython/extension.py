""" Pylada plugin for IPython. """
__all__ = ['load_ipython_extension']

__pylada_is_loaded__ = False
""" Whether the Pylada plugin has already been loaded or not. """
def load_ipython_extension(ip):
  """Load the extension in IPython."""
  global __pylada_is_loaded__
  if not __pylada_is_loaded__:
    from types import ModuleType
    import pylada
    __pylada_is_loaded__ = True
    pylada.interactive = ModuleType('interactive')
    pylada.interactive.jobfolder = None
    pylada.interactive.jobfolder_path = None
    pylada.is_interactive = True

def unload_ipython_extension(ip):
  """ Unloads Pylada IPython extension. """
  ip.user_ns.pop('collect', None)
  ip.user_ns.pop('jobparams', None)
