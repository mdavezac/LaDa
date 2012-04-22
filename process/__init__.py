from ..error import root, out_of_range
class Fail(root):
  """ Process failed to run successfully. """
  pass

def which(program):
  """ Gets location of program by mimicking bash which command. """
  from os import environ
  from os.path import split, expanduser, expandvars, join
  from ..misc import RelativePath
  from ..error import IOError

  def is_exe(path):
    from os import access, X_OK
    from os.path import isfile
    return isfile(path) and access(path, X_OK)

  exprog = expanduser(expandvars(program))
  fpath, fname = split(exprog)
  if fpath:
    if is_exe(exprog): return RelativePath(exprog).path
  else:
    for dir in environ["PATH"].split(':'):
      if is_exe(join(dir, exprog)): return RelativePath(join(dir, exprog)).path

  raise IOError('Could not find executable {0}.'.format(program))
