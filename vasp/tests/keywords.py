def test_bool():
  from lada.vasp import Vasp
  a = Vasp()

  assert a._input['addgrid'].keyword == 'addgrid'
  assert a._input['addgrid'].output_map() is None
  assert a.addgrid is None
  a.addgrid = False
  assert a.addgrid is False
  assert 'addgrid' in a._input['addgrid'].output_map()
  assert a._input['addgrid'].output_map()['addgrid'] == '.FALSE.'
  a.addgrid = True
  assert a.addgrid is True
  assert 'addgrid' in a._input['addgrid'].output_map()
  assert a._input['addgrid'].output_map()['addgrid'] == '.TRUE.'
  a.addgrid = None
  assert a._input['addgrid'].keyword == 'addgrid'
  assert a._input['addgrid'].output_map() is None
  a.addgrid = 0
  assert a.addgrid is False

def test_choice():
  from lada.vasp import Vasp
  from lada.error import ValueError
  a = Vasp()

  assert a.ispin is None
  assert a._input['ispin'].keyword == 'ispin'
  assert a._input['ispin'].output_map() is None
  a.ispin = 1
  assert a.ispin == 1
  assert 'ispin' in a._input['ispin'].output_map()
  assert a._input['ispin'].output_map()['ispin'] == '1'
  a.ispin = 2
  assert a.ispin == 2
  assert 'ispin' in a._input['ispin'].output_map()
  assert a._input['ispin'].output_map()['ispin'] == '2'
  a.ispin = None
  assert a.ispin is None
  assert a._input['ispin'].keyword == 'ispin'
  assert a._input['ispin'].output_map() is None

  try: a.ispin = 5
  except: pass
  else: raise RuntimeError()

  a.ispin = '1'
  assert a.ispin == 1
  a.ispin = '2'
  assert a.ispin == 2

  try: a.ispin = '3'
  except: pass
  else: raise RuntimeError()

def test_alias():
  from lada.vasp import Vasp
  from lada.error import ValueError
  a = Vasp()

  assert a.lmaxmix is None
  assert a._input['lmaxmix'].keyword == 'lmaxmix'
  assert a._input['lmaxmix'].output_map() is None
  for i, channel in enumerate(['s', 'p', 'd', 'f', 'e']):
    a.lmaxmix = channel
    assert a.lmaxmix == channel
    assert 'lmaxmix' in a._input['lmaxmix'].output_map()
    assert a._input['lmaxmix'].output_map()['lmaxmix'] == str(i+1)
    a.lmaxmix = i+1
    assert a.lmaxmix == channel
    a.lmaxmix = str(i+1) 
    assert a.lmaxmix == channel

  try: a.lmaxmix = 'a'
  except ValueError: pass
  else: raise Exception()

def test_typed():
  from lada.vasp import Vasp
  from lada.error import ValueError
  a = Vasp()

  assert a.nbands is None
  assert a._input['nbands'].keyword == 'nbands'
  assert a._input['nbands'].output_map() is None

  a.nbands = 50
  assert a.nbands == 50
  assert 'nbands' in a._input['nbands'].output_map()
  assert a._input['nbands'].output_map()['nbands'] == str(a.nbands)
  a.nbands = '51'
  assert a.nbands == 51
  assert 'nbands' in a._input['nbands'].output_map()
  assert a._input['nbands'].output_map()['nbands'] == str(a.nbands)
  a.nbands = None
  assert a.nbands is None
  assert a._input['nbands'].output_map() is None

  try: a.nbands = 'a'
  except ValueError: pass
  else: raise Exception()

  assert a.smearings is None
  assert a._input['smearings'].keyword == 'smearings'
  assert a._input['smearings'].output_map() is None
  a.smearings = [1.5, 1.0, 0.5]
  assert len(a.smearings) == 3
  assert all(abs(i-v) < 1e-8 for i, v in zip(a.smearings, [1.5, 1.0, 0.5]))
  assert 'smearings' in a._input['smearings'].output_map()
  assert all(abs(float(i)-v) < 1e-8 for i, v in zip(a._input['smearings'].output_map()['smearings'].split(), [1.5, 1.0, 0.5]))
  a.smearings = ['1.2', '0.2']
  assert len(a.smearings) == 2
  assert all(abs(i-v) < 1e-8 for i, v in zip(a.smearings, [1.2, 0.2]))
  assert 'smearings' in a._input['smearings'].output_map()
  assert all(abs(float(i)-v) < 1e-8 for i, v in zip(a._input['smearings'].output_map()['smearings'].split(), [1.2, 0.2]))
  a.smearings = '1.3 0.3'
  assert len(a.smearings) == 2
  assert all(abs(i-v) < 1e-8 for i, v in zip(a.smearings, [1.3, 0.3]))
  assert 'smearings' in a._input['smearings'].output_map()
  assert all(abs(float(i)-v) < 1e-8 for i, v in zip(a._input['smearings'].output_map()['smearings'].split(), [1.3, 0.3]))
  a.smearings = '1.3, 0.3'
  assert len(a.smearings) == 2
  assert all(abs(i-v) < 1e-8 for i, v in zip(a.smearings, [1.3, 0.3]))
  a.smearings = '1.3; 0.3'
  assert len(a.smearings) == 2
  assert all(abs(i-v) < 1e-8 for i, v in zip(a.smearings, [1.3, 0.3]))
  a.smearings = None
  assert a.smearings is None
  assert a._input['smearings'].output_map() is None
  
  try: a.smearings = 5.5
  except ValueError: pass
  else: raise Exception()
  try: a.smearings = [5.5, 'a']
  except ValueError: pass
  else: raise Exception()

if __name__ == '__main__':
  test_bool()
  test_choice()
  test_alias()
  test_typed()
