def raw_input(*args): return 'y'

def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os import makedirs
  from os.path import exists, join
  from IPython.core.interactiveshell import InteractiveShell
  from lada.jobs import JobDict
  import lada
  from dummy import functional
  import __builtin__ 
  __builtin__.raw_input = raw_input

  self = InteractiveShell.instance()

  root = JobDict()
  for type, trial, size in [('this', 0, 10), ('this', 1, 15), ('that', 2, 20), ('that', 1, 20)]:
    job = root / type / str(trial)
    job.functional = functional
    job.params['indiv'] = size
    if type == 'that': job.params['value'] = True

  directory =  mkdtemp() # '/tmp/test' #
  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    self.user_ns['jobdict'] = root
    self.magic("explore jobdict")
    jobdict = lada.interactive.jobdict
    assert 'this/0' in jobdict and 'this/1' in jobdict and 'that/2' in jobdict and 'that/1'
    assert '0' in jobdict['this'] and '1' in jobdict['this']
    assert '1' in jobdict['that'] and '2' in jobdict['that']
    assert 'other' not in jobdict
    for job in jobdict.values():
      assert repr(job.functional) == repr(functional)
    assert getattr(jobdict['this/0'], 'indiv', 0) == 10
    assert getattr(jobdict['this/1'], 'indiv', 0) == 15
    assert getattr(jobdict['that/1'], 'indiv', 0) == 20
    assert getattr(jobdict['that/2'], 'indiv', 0) == 20
    assert not hasattr(jobdict['this/0'], 'value')
    assert not hasattr(jobdict['this/1'], 'value')
    assert getattr(jobdict['that/1'], 'value', False)
    assert getattr(jobdict['that/2'], 'value', False)
    assert lada.interactive.jobdict_path is None
    assert 'jobparams' in self.user_ns
    assert jobdict is self.user_ns['jobparams'].jobdict

    self.magic("savejobs {0}/dict".format(directory))
    lada.interactive.jobdict = None
    lada.interactive.jobdict_path = None
    self.magic("explore {0}/dict".format(directory))
    jobdict = lada.interactive.jobdict
    assert 'this/0' in jobdict and 'this/1' in jobdict and 'that/2' in jobdict and 'that/1'
    assert '0' in jobdict['this'] and '1' in jobdict['this']
    assert '1' in jobdict['that'] and '2' in jobdict['that']
    assert 'other' not in jobdict
    for job in jobdict.values():
      assert repr(job.functional) == repr(functional)
    assert getattr(jobdict['this/0'], 'indiv', 0) == 10
    assert getattr(jobdict['this/1'], 'indiv', 0) == 15
    assert getattr(jobdict['that/1'], 'indiv', 0) == 20
    assert getattr(jobdict['that/2'], 'indiv', 0) == 20
    assert not hasattr(jobdict['this/0'], 'value')
    assert not hasattr(jobdict['this/1'], 'value')
    assert getattr(jobdict['that/1'], 'value', False)
    assert getattr(jobdict['that/2'], 'value', False)
    assert lada.interactive.jobdict_path is not None
    assert 'jobparams' in self.user_ns
    assert jobdict is self.user_ns['jobparams'].jobdict
    assert jobdict is self.user_ns['collect'].jobdict

    for name, job in root.iteritems():
      if name == 'this/1': continue
      job.compute(outdir=join(directory, name))

    self.magic("explore results".format(directory))
    assert set(['/this/0/', '/that/1/', '/that/2/']) == set(self.user_ns['collect'].iterkeys())
    self.magic("explore errors".format(directory))
    assert set(['/this/1/']) == set(self.user_ns['collect'].iterkeys())
    
  finally: 
    if directory != '/tmp/test': rmtree(directory)
    pass
  


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
