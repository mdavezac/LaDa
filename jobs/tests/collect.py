def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os import makedirs
  from os.path import exists, join
  from pickle import dump
  from lada.jobs import JobDict, MassExtract
  from dummy import functional

  root = JobDict()
  for type, trial, size in [('this', 0, 10), ('this', 1, 15), ('that', 2, 20), ('that', 1, 20)]:
    job = root / type / str(trial)
    job.functional = functional
    job.params['indiv'] = size
    if type == 'that': job.params['value'] = True
  job = root / 'this' /  '0' / 'another'
  job.functional = functional
  job.params['indiv'] = 25
  job.params['value'] = 5
  job.params['another'] = 6

  directory =  mkdtemp() # '/tmp/test' # 
  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    for name, job in root.iteritems():
      result = job.compute(outdir=join(directory, name))
      assert result.success
      assert {'this/0': 10, 'this/1': 15, 'that/1': 20, \
              'that/2': 20, 'this/0/another': 25 }[name] == result.indiv
    with open(join(directory, 'dict'), 'w') as file: dump(root, file)
    collect = MassExtract(path=join(directory, 'dict'))

    for i, (name, value) in enumerate(collect.functional.iteritems()):
      assert repr(value) == repr(functional)
    assert i == 4
    for i, (name, value) in enumerate(collect.indiv.iteritems()):
      assert {'/this/0/': 10, '/this/1/': 15, '/that/1/': 20, \
              '/that/2/': 20, '/this/0/another/': 25 }[name] == value,\
             Exception((name, value))
    assert i == 4

    try: collect.that
    except: pass
    else: raise RuntimeError()

    collect.unix_re = True
    def check_items(regex, keys, d):
      i = 0
      for name in d[regex].iterkeys():
        i += 1
        assert name in keys, KeyError((regex, name))
      assert i == len(keys), RuntimeError(regex)
    check_items('*/1', set(['/this/1/', '/that/1/']), collect)
    check_items('this', set(['/this/1/', '/this/0/', '/this/0/another/']), collect)
    check_items('that/2/', set(['/that/2/']),collect)
    job = root / 'this' /  '0' / 'another'
    check_items('*/*/another', ['/this/0/another/'], collect)
    check_items('*/*/another/', ['/this/0/another/'], collect)
    check_items('*/another', [], collect)
    check_items('../that', ['/that/2/', '/that/1/'], collect['this'])
    check_items('../that', [], collect['this/0'])
    check_items('../*', ['/this/0/', '/this/1/', '/this/0/another/'], collect['this/0'])
    check_items('../*', ['/this/0/', '/this/1/', '/this/0/another/'], collect['this/*'])
 
    collect.unix_re = False
    check_items('.*/1', set(['/this/1/', '/that/1/']), collect)
    check_items('this', set(['/this/1/', '/this/0/', '/this/0/another/']), collect)
    check_items('that/2/', set(['/that/2/']), collect)
    check_items('.*/.*/another', ['/this/0/another/'], collect)
    check_items('.*/another', ['/this/0/another/'], collect)
 
    collect.unix_re = True
    collect.naked_end = True
    for i, (key, value) in enumerate(collect['*/1'].indiv.iteritems()):
      assert {'/this/0/': 10, '/this/1/': 15, '/that/1/': 20,
              '/that/2/': 20, '/this/0/another/': 25}[key] == value
    assert i == 1

  finally:
    if directory != '/tmp/test': rmtree(directory)


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
