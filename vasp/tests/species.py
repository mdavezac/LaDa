def test_U():
  """ Test U translation. """
  from pylada.vasp.specie import U
  a = U("liechtenstein", 's', -1e0, 1e0)
  assert a['type'] == 1 and a['l'] == 0 and abs(a['U'] + 1e0) < 1e-8 \
         and abs(a['J'] - 1e0) < 1e-8 and a['func'] == 'U'
  a = U("dudarev", 'p', -0.5, 5e0)
  assert a['type'] == 2 and a['l'] == 1 and abs(a['U'] + 0.5) < 1e-8 \
         and abs(a['J'] - 5e0) < 1e-8 and a['func'] == 'U'
  a = U("dudarev", 'd')
  assert a['type'] == 2 and a['l'] == 2 and abs(a['U']) < 1e-8 \
         and abs(a['J']) < 1e-8 and a['func'] == 'U'
  a = U("dudarev", 'f')
  assert a['type'] == 2 and a['l'] == 3 and abs(a['U']) < 1e-8 \
         and abs(a['J']) < 1e-8 and a['func'] == 'U' 
  try: U('shit')
  except: pass
  else: raise RuntimeError()
  try: U(l=4)
  except: pass
  else: raise RuntimeError()
  try: U(type=0)
  except: pass
  else: raise RuntimeError()

def test_nlep():
  """ Test nlep translation. """
  from pylada.vasp.specie import nlep
  import pylada
  pylada.vasp_has_nlep = True
  a = nlep("liechtenstein", 's', -1e0)
  assert a['type'] == 1 and a['l'] == 0 and abs(a['U0'] + 1e0) < 1e-8 \
         and 'U1' not in a and a['func'] == 'nlep'
  a = nlep("dudarev", 'p', -0.5)
  assert a['type'] == 2 and a['l'] == 1 and abs(a['U0'] + 0.5) < 1e-8 \
         and 'U1' not in a and a['func'] == 'nlep'
  a = nlep("dudarev", 'd')
  assert a['type'] == 2 and a['l'] == 2 and abs(a['U0']) < 1e-8 \
         and 'U1' not in a and a['func'] == 'nlep'
  a = nlep("dudarev", 'f')
  assert a['type'] == 2 and a['l'] == 3 and abs(a['U0']) < 1e-8 \
         and 'U1' not in a and a['func'] == 'nlep'
  try: nlep('shit')
  except: pass
  else: raise RuntimeError()
  try: nlep(l=4)
  except: pass
  else: raise RuntimeError()
  try: nlep(type=0)
  except: pass
  else: raise RuntimeError()

def test_enlep():
  """ Test enlep translation. """
  from pylada.vasp.specie import nlep
  a = nlep("liechtenstein", 's', -1e0, -5e0)
  assert a['type'] == 1 and a['l'] == 0 and abs(a['U0'] + 1e0) < 1e-8 \
         and abs(a['U1'] + 5e0) < 1e-8 and a['func'] == 'enlep'

def test_specie(directory):
  from os.path import join
  from pickle import loads, dumps
  from quantities import eV
  from pylada.vasp.specie import Specie

  pseudos = [('Rh', 229.0, 9.), ('O', 400., 6.), ('Si', 245.345, 4.), ('Zn', 276.727, 12.)]
  for name, enmax, valence in pseudos:
    path = join(join(directory, "pseudos"), name)
    specie = Specie(path)
    specie.potcar_exists()
    assert abs(specie.enmax - enmax*eV) < 1e-8 and abs(specie.valence-valence) < 1e8
    with open(join(path, "POTCAR"), "r") as file:
      specie.read_potcar() == file.read()
    assert repr(specie) == repr(loads(dumps(specie)))


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 2: path.extend(argv[2:])
  
  test_U()
  test_nlep()
  test_enlep()
  if len(argv) > 1: test_specie(argv[1])

