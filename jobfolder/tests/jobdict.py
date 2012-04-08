def test():
  from random import randint
  from tempfile import mkdtemp
  from shutil import rmtree
  from os import makedirs
  from os.path import exists, join
  from pickle import loads, dumps
  from lada.jobs import JobFolder
  from dummy import functional

  sizes = [10, 15, 20, 25]
  root = JobFolder()
  for type, trial, size in [('this', 0, 10), ('this', 1, 15), ('that', 2, 20), ('that', 1, 20)]:
    job = root / type / str(trial)
    job.functional = functional
    job.params['indiv'] = size
    if type == 'that': job.params['value'] = True
      
  assert 'this/0' in root and 'this/1' in root and 'that/2' in root and 'that/1'
  assert '0' in root['this'] and '1' in root['this']
  assert '1' in root['that'] and '2' in root['that']
  assert 'other' not in root
  for job in root.values():
    assert repr(job.functional) == repr(functional)
  assert getattr(root['this/0'], 'indiv', 0) == 10
  assert getattr(root['this/1'], 'indiv', 0) == 15
  assert getattr(root['that/1'], 'indiv', 0) == 20
  assert getattr(root['that/2'], 'indiv', 0) == 20
  assert not hasattr(root['this/0'], 'value')
  assert not hasattr(root['this/1'], 'value')
  assert getattr(root['that/1'], 'value', False)
  assert getattr(root['that/2'], 'value', False)

  for key, test in zip(root, ['that/1', 'that/2', 'this/0', 'this/1']): 
    assert key == test
  for key, test in zip(root['that/1'].root, ['that/1', 'that/2', 'this/0', 'this/1']): 
    assert key == test
  for key, test in zip(root['that'], ['1', '2']): assert key == test
  for key, test in zip(root['this'], ['0', '1']): assert key == test
  del root['that/2']
  assert 'that/2' not in root
  
  directory = mkdtemp() # '/tmp/test'
  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    for name, job in root.iteritems():
      result = job.compute(outdir=join(directory, name))
      assert result.success
      assert exists(join(directory, name))
      assert exists(join(join(directory, name), 'OUTCAR'))
      if name == 'that/1': assert result.indiv == 20
      elif name == 'that/2': assert result.indiv == 15
      elif name == 'this/1': assert result.indiv == 15
      elif name == 'this/0': assert result.indiv == 10
  finally:
    if directory != '/tmp/test': rmtree(directory)

  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    for name, job in loads(dumps(root)).iteritems():
      result = job.compute(outdir=join(directory, name))
      assert result.success
      assert exists(join(directory, name))
      assert exists(join(join(directory, name), 'OUTCAR'))
      if name == 'that/1': assert result.indiv == 20
      elif name == 'that/2': assert result.indiv == 15
      elif name == 'this/1': assert result.indiv == 15
      elif name == 'this/0': assert result.indiv == 10
  finally:
    if directory != '/tmp/test': rmtree(directory)


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
