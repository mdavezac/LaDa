def test():
  from pickle import loads, dumps
  from lada.vasp.incar._params import Algo
  import lada
  lada.is_vasp_4 = True

  # default.
  a = Algo()
  assert a.value == "Fast"
  # wrong argument.
  try: a.value = 0
  except: pass
  else: raise RuntimeError()
  try: a.value = "WTF"
  except: pass
  else: raise RuntimeError()

  # possible inputs and some.
  d = {
      'Very_Fast': ['very fast', 'VERY-fAst', 'very_FAST', 'v'],
      'VeryFast': ['very fast', 'VERY-fAst', 'very_FAST', 'v'],
      'Fast': ['fast', 'f'],
      'Normal': ['normal', 'n'],
      'Damped': ['damped', 'd'],
      'Diag': ['diag'],
      'All': ['all', 'a'],
      'Nothing': ['nothing'],
      'chi': ['chi'],
      'GW': ['gw'],
      'GW0': ['gw0'],
      'scGW': ['scgw'],
      'scGW0': ['scgw0'],
      'Conjugate': ['conjugate', 'c'],
      'Subrot': ['subrot', 's'],
      'Eigenval': ['eigenval', 'e']
  }
  vasp5 = 'Subrot', 'chi', 'GW', 'GW0', 'scGW', 'scGW0', 'Conjugate', 'Eigenval', 'Exact', 'Nothing'
  for isvasp4 in [True, False]:
    lada.is_vasp_4 = isvasp4
    for key, items in d.iteritems():
      for value in items:
        if key in vasp5 and isvasp4:
          try: a.value = value
          except: pass 
          else: raise RuntimeError((value, key))
        elif key == 'VeryFast' and isvasp4:
          a.value = value
          assert a.value == 'Very_Fast'
          assert loads(dumps(a)).incar_string() == "ALGO = {0}".format('Very_Fast')
          assert repr(a) == "Algo('Very_Fast')"
        elif key == 'Very_Fast' and not isvasp4:
          a.value = value
          assert a.value == 'VeryFast'
          assert loads(dumps(a)).incar_string() == "ALGO = {0}".format('VeryFast')
          assert repr(a) == "Algo('VeryFast')"
        else:
          a.value = value
          assert a.value == key
          assert loads(dumps(a)).incar_string() == "ALGO = {0}".format(key)
          assert repr(a) == "Algo({0!r})".format(key)


if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

