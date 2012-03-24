# import _doc
# __doc__ = _doc.__doc__
__docformat__ = "restructuredtext en"

__lada_is_loaded__ = False
""" Whether the LaDa plugin has already been loaded or not. """
def load_ipython_extension(ip):
  """Load the extension in IPython."""
  global __lada_is_loaded__
  if not __lada_is_loaded__:
    from types import ModuleType, MethodType
    from .savejobs import savejobs
    from .explore import explore, completer as explore_completer
    from .goto import goto, completer as goto_completer
    from .listjobs import listjobs
    from .showme import showme, completer as showme_completer
    import lada
    __lada_is_loaded__ = True
    lada.interactive = ModuleType('interactive')
    lada.interactive.jobdict = None
    lada.interactive.jobdict_path = None
    lada.is_interactive = True
    ip.define_magic('savejobs', savejobs)
    ip.define_magic('explore', explore)
    ip.define_magic('goto', goto)
    ip.define_magic('listjobs', listjobs)
    ip.define_magic('showme', showme)
    ip.set_hook('complete_command', explore_completer, str_key = '%explore')
    ip.set_hook('complete_command', goto_completer, str_key = '%goto')
    ip.set_hook('complete_command', showme_completer, str_key = '%showme')

def unload_ipython_extension(ip):
  """ Unloads LaDa IPython extension. """
  ip.user_ns.pop('collect', None)
  ip.user_ns.pop('jobparams', None)

# def cancel_completer(self, info):
#   return qstat(self, info.symbol).fields(-1)[1:]

# def cancel_jobs(self, arg):
#   """ Cancel jobs which grep for whatever is in arg.
#   
#       For instance, the following cancels all jobs with "anti-ferro" in their
#       name.

#       >>> %cancel_jobs "anti-ferro"
#   """
#   from lada import lada_with_slurm
#   from .qstat import qstat
#   arg = str(arg[1:-1])
#   if len(arg) != 0: 
#     result = qstat(self, arg)
#     for u, name in zip(result.fields(0), result.fields(-1)):
#       print "cancelling %s." % (name)
#     message = "Are you sure you want to cancel the jobs listed above? [y/n] "
#   else: message = "Cancel all jobs? [y/n] "
#   a = ''
#   while a not in ['n', 'y']: a = raw_input(message)
#   if a == 'n': return
#   
#   cmd = "scancel " if lada_with_slurm  else  "qdel "
#   result = qstat(self, arg)
#   for u, name in zip(result.fields(0), result.fields(-1)): self.api.system(cmd + str(u))

