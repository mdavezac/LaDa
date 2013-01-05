def raw_input(*args): return 'y'

def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os import makedirs, getcwd, chdir
  from os.path import exists, join
  from IPython.core.interactiveshell import InteractiveShell
  from pylada.jobfolder import JobFolder
  from pylada import interactive
  from dummy import functional
  import __builtin__ 
  try: 
    saveri = __builtin__.raw_input
    __builtin__.raw_input = raw_input
 
    self = InteractiveShell.instance()
 
    root = JobFolder()
    for type, trial, size in [('this', 0, 10), ('this', 1, 15), ('that', 2, 20), ('that', 1, 20)]:
      job = root / type / str(trial)
      job.functional = functional
      job.params['indiv'] = size
      if type == 'that': job.params['value'] = True
 
    origdir = getcwd()
    directory = mkdtemp()
    if exists(directory) and directory == '/tmp/test': rmtree(directory)
    if not exists(directory): makedirs(directory)
    try: 
      self.user_ns['jobfolder'] = root
      self.magic("explore jobfolder")
      self.magic("savefolders {0}/dict".format(directory))
      for name, job in root.iteritems():
        result = job.compute(outdir=join(directory, name))
        assert result.success
        assert {'this/0': 10, 'this/1': 15, 'that/1': 20, \
                'that/2': 20, 'this/0/another': 25 }[name] == result.indiv
 
      self.magic("explore {0}/dict".format(directory))
      self.magic("goto this/0")
      assert getcwd() == '{0}/this/0'.format(directory)
      assert interactive.jobfolder.name == '/this/0/'
      self.magic("goto ../1")
      assert getcwd() == '{0}/this/1'.format(directory)
      assert interactive.jobfolder.name == '/this/1/'
      self.magic("goto /that")
      assert getcwd() == '{0}/that'.format(directory)
      assert interactive.jobfolder.name == '/that/'
      self.magic("goto 2")
      assert getcwd() == '{0}/that/2'.format(directory)
      assert interactive.jobfolder.name == '/that/2/'
      self.magic("goto /")
      self.magic("goto next")
      assert getcwd() == '{0}/that/1'.format(directory)
      assert interactive.jobfolder.name == '/that/1/'
      self.magic("goto next")
      assert getcwd() == '{0}/that/2'.format(directory)
      assert interactive.jobfolder.name == '/that/2/'
      self.magic("goto previous")
      assert getcwd() == '{0}/that/1'.format(directory)
      assert interactive.jobfolder.name == '/that/1/'
      self.magic("goto next")
      assert getcwd() == '{0}/this/0'.format(directory)
      assert interactive.jobfolder.name == '/this/0/'
      self.magic("goto next")
      assert getcwd() == '{0}/this/1'.format(directory)
      assert interactive.jobfolder.name == '/this/1/'
      self.magic("goto next") # no further jobs
      assert getcwd() == '{0}/this/1'.format(directory)
      assert interactive.jobfolder.name == '/this/1/'
      self.magic("goto /") # go back to root to avoid errors
      
    finally: 
      chdir(origdir)
      try: 
        if directory != '/tmp/test': rmtree(directory)
      except: pass
  finally: __builtin__.raw_input = saveri


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
