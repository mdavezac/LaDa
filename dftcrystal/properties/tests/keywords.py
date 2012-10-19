def test_rdfmwf():
  from collections import namedtuple
  from shutil import rmtree
  from tempfile import mkdtemp
  from os.path import join, exists
  from os import remove
  from lada.dftcrystal.properties import Properties
  from lada.error import IOError

  Extract = namedtuple('Extract', ['directory'])
  try: 
    indir = mkdtemp()
    outdir = mkdtemp()

    extract = Extract(indir)
    a = Properties(input=extract)
    assert a.rdfmwf is None
    assert len(a.output_map(filework=False)) == 0
    assert len(a.print_input(filework=False)) == 0
    try: a.output_map(filework=True, workdir=outdir)
    except IOError: pass
    else: raise Exception()

    a.rdfmwf = True
    assert len(a.output_map(filework=False)) == 1
    assert a.output_map(filework=False)['rdfmwf'] == True
    try: a.output_map(filework=True, workdir=outdir)
    except IOError: pass
    else: raise Exception()

    a.rdfmwf = False
    assert len(a.output_map(filework=False)) == 0
    try: a.output_map(filework=True, workdir=outdir)
    except IOError: pass
    else: raise Exception()

    # now with .f98 file
    with open(join(indir, 'crystal.f98'), 'w') as file: file.write('hello')
    a.rdfmwf = None
    assert len(a.output_map(filework=False)) == 1
    assert a.output_map(filework=False)['rdfmwf'] == True
    assert not exists(join(outdir, 'crystal.f98'))
    # should fail
    a.rdfmwf = False
    try: a.output_map(filework=True, workdir=outdir)
    except IOError: pass
    else: raise Exception()
    # now check writing
    a.rdfmwf = None
    assert not exists(join(outdir, 'crystal.f98'))
    assert a.output_map(filework=True, workdir=outdir)['rdfmwf'] == True
    assert exists(join(outdir, 'crystal.f98'))
    remove(join(outdir, 'crystal.f98'))
    a.rdfmwf = True
    assert not exists(join(outdir, 'crystal.f98'))
    assert a.output_map(filework=True, workdir=outdir)['rdfmwf'] == True
    assert exists(join(outdir, 'crystal.f98'))
    remove(join(indir, 'crystal.f98'))
    assert exists(join(outdir, 'crystal.f98'))
    assert a.output_map(filework=True, workdir=outdir, outdir=outdir)['rdfmwf'] == True
    assert exists(join(outdir, 'crystal.f98'))
    remove(join(outdir, 'crystal.f98'))
    
    # now for .f9
    with open(join(indir, 'crystal.f9'), 'w') as file: file.write('hello')
    a.rdfmwf = None
    assert len(a.output_map(filework=False)) == 0
    assert not exists(join(outdir, 'crystal.f9'))
    # should fail
    a.rdfmwf = True
    try: a.output_map(filework=True, workdir=outdir)
    except IOError: pass
    else: raise Exception()
    # now check writing
    a.rdfmwf = None
    assert not exists(join(outdir, 'crystal.f9'))
    assert len(a.output_map(filework=True, workdir=outdir)) == 0
    assert exists(join(outdir, 'crystal.f9'))
    remove(join(outdir, 'crystal.f9'))
    a.rdfmwf = False
    assert not exists(join(outdir, 'crystal.f9'))
    assert len(a.output_map(filework=True, workdir=outdir)) == 0
    assert exists(join(outdir, 'crystal.f9'))
    remove(join(indir, 'crystal.f9'))
    assert exists(join(outdir, 'crystal.f9'))
    assert len(a.output_map(filework=True, workdir=outdir, outdir=outdir)) == 0
    assert exists(join(outdir, 'crystal.f9'))
    remove(join(outdir, 'crystal.f9'))

  finally:
    try: rmtree(indir)
    except: pass
    try: rmtree(outdir)
    except: pass

def test_nosymada():
  from lada.dftcrystal.properties import Properties
  a = Properties()
  assert a.nosymada is None
  assert len(a.output_map()) == 0
  a.nosymada = True
  assert a.print_input() == "NOSYMADA\nEND\n"

def test_newk():
  from pickle import loads, dumps
  from numpy import array, all
  from lada.dftcrystal.properties import Properties
  from lada.dftcrystal.properties.keywords import NewK
  from lada.dftcrystal.input import SetPrint
  a = Properties()
  assert len(a.output_map()) == 0
  a.newk.mp = 5
  assert len(a.output_map()) != 0
  assert all(array(a.output_map()['newk'].split(), dtype='int64') == [5, 5, 0, 0])
  a.newk.printing[66] = -5
  assert all(array(a.output_map()['newk'].split(), dtype='int64') == [5, 5, 0, 1, 66, -5])

  b = a._input['newk']
  assert eval(repr(b), {'NewK': NewK}).output_map()['newk'] == a.output_map()['newk']
  assert loads(dumps(b)).output_map()['newk'] == a.output_map()['newk']
  assert b.__ui_repr__({})['newk'] == repr(b)
  assert b.__ui_repr__({}, defaults=NewK())['newk'] == repr(b)
  del b.printing
  assert len(b.printing) == 0
  assert type(b.printing) is SetPrint
  assert b.__ui_repr__({})['newk'] == repr(b)
  assert b.__ui_repr__({}, defaults=NewK())['newk'] == '5, None'

  assert not b.recompute_fermi
  a.newk.recompute_fermi = True
  assert b.recompute_fermi
  assert all(array(a.output_map()['newk'].split(), dtype='int64') == [5, 5, 1, 0])
  assert eval(repr(b), {'NewK': NewK}).output_map()['newk'] == a.output_map()['newk']
  assert loads(dumps(b)).output_map()['newk'] == a.output_map()['newk']
  assert b.__ui_repr__({})['newk'] == repr(b)
  assert b.__ui_repr__({}, defaults=NewK())['newk'] == repr(b)

def test_band():
  from pickle import loads, dumps
  from collections import namedtuple
  from random import randint
  from numpy import all, abs, array
  from lada.dftcrystal.properties import Properties
  from lada.dftcrystal.properties.keywords import Band

  Shell = namedtuple('Shell', ['charge'])
  Functional = namedtuple('Functional', ['basis'])
  Extract    = namedtuple('Extract', ['functional', 'directory'])
  Atom       = namedtuple('Atom', ['type'])
  class Structure(list):
    def __init__(self, *args, **kwargs):
      super(Structure, self).__init__(*args, **kwargs)
    def eval(self): return self
  shells     = {}
  shells['Al'] = [Shell(randint(0, 10)) for i in xrange(5)]
  shells['Be'] = [Shell(randint(0, 10)) for i in xrange(5)]
  functional = Functional(shells)
  input = Extract(functional, '/dev/null')
  structure = Structure([Atom(u) for u in ['Al']*3 + ['Be'] * 4])
  kwargs = {'input': input, 'structure': structure, 'filework': False}

  nelectrons = sum(u.charge for u in shells['Al']) * 3                         \
               + sum(u.charge for u in shells['Be']) * 4  

  a = Properties()
  assert len(a.output_map()) == 0
  assert len(a.output_map(**kwargs)) == 0
  a.band += [[0, 1, 0], [0, 0, 0]], [[0, 0, 0], [1, 0, 0]]
  assert all(abs(eval(repr(a.band), {'Band': Band})._lines - a.band._lines) < 1e-8)
  assert all(abs(loads(dumps(a.band))._lines - a.band._lines) < 1e-8)
  c = a.output_map(**kwargs)['band'].split()
  assert all(array(c[:7], dtype='float64') == [2, 100000000, 40, 0, nelectrons+6, 1, 0])
  assert all(abs(array(c[7:], dtype='float64') / 1e8 - a.band._lines.flatten()) < 1e-8)
  assert a._input['band'].__ui_repr__({}, name='band')['band'] == repr(a.band)
  c = a._input['band'].__ui_repr__({}, name='band', defaults=Band())['band']
  assert all(abs(array(eval(c), dtype='float64') - a.band._lines) < 1e-8)
  
  a.band.maxband = nelectrons // 2 + 1
  assert a.band.maxband == nelectrons // 2 + 1
  assert all(abs(eval(repr(a.band), {'Band': Band})._lines - a.band._lines) < 1e-8)
  assert eval(repr(a.band), {'Band': Band}).maxband == nelectrons // 2 + 1
  assert all(abs(loads(dumps(a.band))._lines - a.band._lines) < 1e-8)
  assert loads(dumps(a.band)).maxband == nelectrons // 2 + 1
  assert a._input['band'].__ui_repr__({}, name='band')['band'] == repr(a.band)
  assert a._input['band'].__ui_repr__({}, name='band', defaults=Band())['band'] == repr(a.band)
  
  # check title
  assert len(a.output_map(properties=a, **kwargs)['band'].splitlines()[0]) == 0
  structure.name = 'hello'
  assert a.output_map(properties=a, **kwargs)['band'].splitlines()[0] == 'hello'
  a.band.title = 'goodbye'
  assert a.output_map(properties=a, **kwargs)['band'].splitlines()[0] == 'goodbye'

if __name__ == '__main__': 
  test_rdfmwf()
  test_nosymada()
  test_newk()
  test_band()
