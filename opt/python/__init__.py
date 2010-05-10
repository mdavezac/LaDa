""" Miscellaneous """
import _opt 
__load_vasp_in_global_namespace__ = _opt.__load_vasp_in_global_namespace__
__load_pescan_in_global_namespace__ = _opt.__load_pescan_in_global_namespace__
cReals = _opt.cReals
_RedirectFortran = _opt._RedirectFortran
streams = _opt._RedirectFortran.fortran

class _RedirectPy:
  """ Redirects C/C++ input, output, error. """
  def __init__(self, unit, filename, append):
    self.unit = unit
    self.filename = filename
    self.append = append
  def __enter__(self):
    import sys
    if self.unit == streams.input:    self.old = sys.stdin
    elif self.unit == streams.error:  self.old = sys.stderr
    elif self.unit == streams.output: self.old = sys.stdout
    else: raise RuntimeError("Unknown redirection unit.")
    self.file = open(self.filename if len(self.filename) else "/dev/null", "a" if self.append else "w")
    if self.unit == streams.input:    sys.stdin  = self.file
    elif self.unit == streams.error:  sys.stderr = self.file
    elif self.unit == streams.output: sys.stdout = self.file
    else: raise RuntimeError("Unknown redirection unit.")
    return self
  def __exit__(self, *wargs):
    import sys 
    if self.unit == streams.input:    sys.stdin  = self.old 
    elif self.unit == streams.error:  sys.stderr = self.old 
    elif self.unit == streams.output: sys.stdout = self.old 
    else: raise RuntimeError("Unknown redirection unit.")
    self.file.close()
      
class _RedirectAll:
  """ Redirects C/C++ input, output, error. """
  def __init__(self, unit, filename, append):
    self.unit = unit
    self.filename = filename
    if self.filename == None: self.filename = ""
    if self.filename == "/dev/null": self.filename = ""
    self.append = append
  def __enter__(self):
    import os
    import sys
    from boost.mpi import world as comm
    if self.unit == streams.input:    self.old = os.dup(sys.__stdin__.fileno())
    elif self.unit == streams.output: self.old = os.dup(sys.__stdout__.fileno())
    elif self.unit == streams.error:  self.old = os.dup(sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")
    self.filefd = os.open(self.filename if len(self.filename) else os.devnull,\
        (os.O_RDWR | os.O_APPEND) if self.append else (os.O_CREAT | os.O_WRONLY | os.O_EXCL))
    if self.append: self.file = os.fdopen(self.filefd, 'a')
    else: self.file = os.fdopen(self.filefd, 'w')
    if self.unit == streams.input:    os.dup2(self.file.fileno(), sys.__stdin__.fileno())
    elif self.unit == streams.output: os.dup2(self.file.fileno(), sys.__stdout__.fileno())
    elif self.unit == streams.error:  os.dup2(self.file.fileno(), sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")
    return self
  def __exit__(self, *wargs):
    import os
    import sys
    try: self.file.close()
    except: print "Am here"
    if self.unit == streams.input:    os.dup2(self.old, sys.__stdin__.fileno()) 
    elif self.unit == streams.output: os.dup2(self.old, sys.__stdout__.fileno())
    elif self.unit == streams.error:  os.dup2(self.old, sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")

def redirect_all(output=None, error=None, input=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      @param output: Filename to which to redirect output. 
      @param error: Filename to which to redirect err. 
      @param input: Filename to which to redirect input. 
      @param append: If true, will append to files. All or nothing.
  """
  from contextlib import nested
  result = []
  for value, unit in [ (output, streams.output), (error, streams.error), (input, streams.input) ]:
    if value == None: continue
    result.append( _RedirectAll(unit=unit, filename=value, append=append) )
  return nested(*result)


def redirect(fout=None, ferr=None, fin=None, pyout=None, pyerr=None, pyin=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      @param fout: Filename to which to redirect fortran output. 
      @param ferr: Filename to which to redirect fortran err. 
      @param fin: Filename to which to redirect fortran input. 
      @param pyout: Filename to which to redirect C/C++ output. 
      @param pyerr: Filename to which to redirect C/C++ err. 
      @param pyin: Filename to which to redirect C/C++ input. 
      @param append: If true, will append to files. All or nothing.
  """
  from contextlib import nested
  result = []
  for value, unit in [ (fout, streams.output), (ferr, streams.error), (fin, streams.input) ]:
    if value == None: continue
    result.append( _RedirectFortran(unit=unit, filename=value, append=append) )
  for value, unit in [ (pyout, streams.output), (pyerr, streams.error), (pyin, streams.input) ]:
    if value == None: continue
    result.append( _RedirectPy(unit=unit, filename=value, append=append) )
  return nested(*result)


