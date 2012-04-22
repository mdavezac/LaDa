""" LaDa plugin for IPython. """
__all__ = ['load_ipython_extension']

__lada_is_loaded__ = False
""" Whether the LaDa plugin has already been loaded or not. """
def load_ipython_extension(ip):
  """Load the extension in IPython."""
  global __lada_is_loaded__
  if not __lada_is_loaded__:
    from types import ModuleType
    import lada
    __lada_is_loaded__ = True
    lada.interactive = ModuleType('interactive')
    lada.interactive.jobfolder = None
    lada.interactive.jobfolder_path = None
    lada.is_interactive = True

def unload_ipython_extension(ip):
  """ Unloads LaDa IPython extension. """
  ip.user_ns.pop('collect', None)
  ip.user_ns.pop('jobparams', None)
