def test():
  from random import randint
  from tempfile import mkdtemp
  from shutil import rmtree
  from os import makedirs
  from os.path import exists, join
  from lada.jobs import JobDict, JobParams
  from dummy import functional

  sizes = [10, 15, 20, 25]
  root = JobDict()
  for type, trial, size in [('this', 0, 10), ('this', 1, 15), ('that', 2, 20), ('that', 1, 20)]:
    job = root / type / str(trial)
    job.functional = functional
    job.params['indiv'] = size
    if type == 'that': job.params['value'] = True

  jobparams = JobParams(root)
  assert len(jobparams.functional.values()) == 4
  for i, (name, value) in enumerate(jobparams.functional.iteritems()):
    assert repr(value) == repr(functional)
  assert i == 3
  for i, (name, value) in enumerate(jobparams.indiv.iteritems()):
    if   name == '/this/0/': assert value == 10
    elif name == '/this/1/': assert value == 15
    elif name == '/that/1/': assert value == 20
    elif name == '/that/2/': assert value == 20
    else: raise RuntimeError()
  assert i == 3

  for i, (name, value) in enumerate(jobparams.value.iteritems()):
    if   name == '/that/1/': assert value == True
    elif name == '/that/2/': assert value == True
    else: raise RuntimeError()
  assert i == 1
  
  try: jobparams.that
  except: pass
  else: raise RuntimeError()

  jobparams.unix_re = True
  def check_items(regex, keys, d):
    i = 0
    for name in d[regex].iterkeys():
      i += 1
      assert name in keys, KeyError((regex, name))
    assert i == len(keys), RuntimeError(regex)
  check_items('*/1', set(['/this/1/', '/that/1/']), jobparams)
  check_items('this', set(['/this/1/', '/this/0/']), jobparams)
  check_items('that/2/', set(['/that/2/']),jobparams)
  job = root / 'this' /  '0' / 'another'
  job.functional = functional
  job.params['indiv'] = 25
  job.params['value'] = 5
  job.params['another'] = 6
  check_items('*/*/another', ['/this/0/another/'], jobparams)
  check_items('*/*/another/', ['/this/0/another/'], jobparams)
  check_items('*/another', [], jobparams)
  check_items('../that', ['/that/2/', '/that/1/'], jobparams['this'])
  check_items('../that', [], jobparams['this/0'])
  check_items('../*', ['/this/0/', '/this/1/', '/this/0/another/'], jobparams['this/0'])
  check_items('../*', ['/this/0/', '/this/1/', '/this/0/another/'], jobparams['this/*'])

  jobparams.unix_re = False
  check_items('.*/1', set(['/this/1/', '/that/1/']), jobparams)
  check_items('this', set(['/this/1/', '/this/0/', '/this/0/another/']), jobparams)
  check_items('that/2/', set(['/that/2/']), jobparams)
  check_items('.*/.*/another', ['/this/0/another/'], jobparams)
  check_items('.*/another', ['/this/0/another/'], jobparams)

  jobparams.unix_re = True
  jobparams.naked_end = True
  assert isinstance(jobparams.another, int)
  for i, (key, value) in enumerate(jobparams['*/1'].indiv.iteritems()):
    assert {'/this/0/': 10, '/this/1/': 15, '/that/1/': 20,
            '/that/2/': 20, '/this/0/another/': 25}[key] == value
  assert i == 1

  jobparams.indiv = 5
  for i, (key, value) in enumerate(jobparams.indiv.iteritems()): assert value == 5
  assert i == 4
  jobparams['*/0'].indiv = 6
  assert len(jobparams.indiv.values()) == 5
  for i, (key, value) in enumerate(jobparams.indiv.iteritems()): 
    assert value == (6 if '0' in key else 5)
  assert i == 4

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
