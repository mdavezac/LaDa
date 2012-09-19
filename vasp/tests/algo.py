def test():
  from pickle import loads, dumps
  from lada.vasp import Vasp
  import lada
  
  lada.is_vasp_4 = True
  a = Vasp()

  # default.
  assert a.algo == 'Fast'
  # wrong argument.
  try: a.algo = 0
  except: pass
  else: raise RuntimeError()
  try: a.algo = "WTF"
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
  dictionary = {'Algo': a._input['algo'].__class__}
  for isvasp4 in [True, False]:
    lada.is_vasp_4 = isvasp4
    for key, items in d.iteritems():
      for value in items:
        if key in vasp5 and isvasp4:
          try: a.algo = value
          except: pass 
          else: raise RuntimeError((value, key))
          continue
        a.algo = value
        o = a._input['algo']
        if key == 'VeryFast' and isvasp4:
          assert a.algo == 'Very_Fast'
          assert o.output_map()['algo'] == 'Very_Fast'
          assert loads(dumps(o)).output_map()["algo"] == 'Very_Fast'
          assert eval(repr(o), dictionary).output_map()["algo"] == 'Very_Fast'
        elif key == 'Very_Fast' and not isvasp4:
          assert a.algo == 'VeryFast'
          assert o.output_map()['algo'] == 'VeryFast'
          assert loads(dumps(o)).output_map()["algo"] == 'VeryFast'
          assert eval(repr(o), dictionary).output_map()["algo"] == 'VeryFast'
        else:
          assert a.algo == key
          assert o.output_map()['algo'] == key
          assert loads(dumps(o)).output_map()["algo"] == key
          assert eval(repr(o), dictionary).output_map()["algo"] == key

  a.algo = None
  assert a.algo is None
  assert o.output_map() is None
  assert loads(dumps(o)).output_map() is None
  assert eval(repr(o), dictionary).output_map() is None


if __name__ == '__main__': 
  test()

 
