# import _doc
# __doc__ = _doc.__doc__
__docformat__ = "restructuredtext en"

__lada_is_loaded__ = False
""" Whether the LaDa plugin has already been loaded or not. """
def load_ipython_extension(ip):
  """Load the extension in IPython."""
  global __lada_is_loaded__
  if not __lada_is_loaded__:
    from types import ModuleType
    from sys import modules
    from .savefolders import savefolders
    from .explore import explore, completer as explore_completer
    from .goto import goto, completer as goto_completer
    from .listfolders import listfolders
    from .showme import showme, completer as showme_completer
    from .launch import launch, completer as launch_completer
    from .export import export, completer as export_completer
    import lada
    __lada_is_loaded__ = True
    lada.interactive = ModuleType('interactive')
    lada.interactive.jobfolder = None
    lada.interactive.jobfolder_path = None
    lada.is_interactive = True
    modules['lada.interactive'] = lada.interactive
    ip.define_magic('savefolders', savefolders)
    ip.define_magic('explore', explore)
    ip.define_magic('goto', goto)
    ip.define_magic('listfolders', listfolders)
    ip.define_magic('showme', showme)
    ip.define_magic('launch', launch)
    ip.define_magic('export', export)
    ip.set_hook('complete_command', explore_completer, str_key = '%explore')
    ip.set_hook('complete_command', goto_completer, str_key = '%goto')
    ip.set_hook('complete_command', showme_completer, str_key = '%showme')
    ip.set_hook('complete_command', launch_completer, str_key = '%launch')
    ip.set_hook('complete_command', export_completer, str_key = '%export')
    if lada.ipython_verbose_representation is not None:
      lada.verbose_representation = lada.ipython_verbose_representation
    if hasattr(lada, 'ipython_qstat'):
      ip.define_magic('qstat', qstat)
      ip.define_magic('qdel', qdel)
      def dummy(*args, **kwargs): return []
      ip.set_hook('complete_command', dummy, str_key = '%qdel')
      ip.set_hook('complete_command', dummy, str_key = '%qstat')

def unload_ipython_extension(ip):
  """ Unloads LaDa IPython extension. """
  ip.user_ns.pop('collect', None)
  ip.user_ns.pop('jobparams', None)

def get_shell(self):
  """ Gets shell despite ipython version issues """
  return getattr(self, 'shell', self)

def jobfolder_file_completer(data):
  """ Returns list of potential job-folder and directories. """
  from os.path import isdir
  from glob import iglob
  from IPython.core.completer import expand_user, compress_user
  from .. import jobfolder_glob
  if len(data) == 0: data = ['']
  relpath, tilde_expand, tilde_val = expand_user(data[-1])
  dirs = [f.replace('\\','/') + "/" for f in iglob(relpath+'*') if isdir(f)]
  dicts = [ f.replace('\\','/') for u in jobfolder_glob for f in iglob(relpath+u)]
  if '.' in data[-1]:
    relpath, a, b = expand_user(data[-1][:data[-1].find('.')])
    dicts.extend([ f.replace('\\','/') for u in jobfolder_glob for f in iglob(relpath+u)])
  dummy = [compress_user(p, tilde_expand, tilde_val) for p in dirs+dicts]
  return [d for d in dummy if d not in data]


def save_n_explore(folder, path):
  """ Save and explore job-folder. 

      For use with ipython interactive terminal only.
  """ 
  from .. import is_interactive
  from ..error import interactive as ierror
  if not is_interactive: raise ierror('Not in interactive session.')

  from IPython.core.interactiveshell import InteractiveShell
  from ..ipython.explore import explore
  from ..ipython.savefolders import savefolders
  import lada
  
  lada.interactive.jobfolder = folder.root
  lada.interactive.jobfolder_path = path
  shell = InteractiveShell.instance()
  savefolders(shell, path)
  explore(shell, '{0}  --file'.format(path))

def qdel_completer(self, info):
  """ Completer for qdel. 

      Too slow. Disabled.
  """
  return self.magic("%qstat {0}".format(info.symbol)).fields(-1)

def qdel(self, arg):
  """ Cancel jobs which grep for whatever is in arg.
  
      For instance, the following cancels all jobs with "anti-ferro" in their
      name. The name is the last column in qstat.

      >>> %qdel "anti-ferro"
  """
  from lada import qdel_exe
  arg = arg.lstrip().rstrip()
  if '--help' in arg.split() or '-h' in arg.split():
   print qdel.__doc__
   return 
  if len(arg) != 0: 
    result = self.qstat(arg)
    if len(result) == 0:
      print 'No jobs in queue'
      return
    for u, name in zip(result.fields(0), result.fields(-1)):
      print "cancelling %s." % (name)
    message = "Are you sure you want to cancel the jobs listed above? [y/n] "
  else: message = "Cancel all jobs? [y/n] "
  a = ''
  while a not in ['n', 'y']: a = raw_input(message)
  if a == 'n': return
  
  result = qstat(self, arg)
  for u, name in zip(result.fields(0), result.fields(-1)): 
    self.shell.system('{0} {1}'.format(qdel_exe, u))

def qstat(self, arg):
  """ SList of user's jobs. 

      The actual print-out and functionality will depend on the user-specified
      function :py:func:`lada.ipython_qstat`. However, in general %qstat should
      work as follows:
 
      >>> %qstat
      [ 'id something something jobname' ]
      
      It returns an SList_ of all the users jobs, with the job-id as the first
      column and the job-name ass the last. The results can be filtered using
      SList_'s grep, or directly as in:

      >>> %qstat hello
      [ 'id something something hellorestofname' ]

      .. _SList: http://ipython.org/ipython-doc/stable/api/generated/IPython.utils.text.html#slist
  """
  from lada import ipython_qstat
  arg = arg.rstrip().lstrip()
  if len(arg) != 0 and '--help' in arg.split() or '-h' in arg.split():
    print qstat.__doc__ + '\n' + ipython_qstat.__doc__
    return
  result = ipython_qstat(self, arg)
  if len(arg) == 0: return result
  return result.grep(arg, field=-1)
