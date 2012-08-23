def test_shrink():
  from pickle import loads, dumps
  from numpy import array, all
  from lada.dftcrystal.electronic import Shrink
  from lada.dftcrystal.molecule import Molecule

  a = Shrink()
  assert a.mp == 1
  assert a.gallat is None
  assert 'shrink' in a.output_map()
  assert a.output_map()['shrink'] == a.raw
  assert a.output_map()['shrink'] == '1 1\n'
  assert repr(a) == 'Shrink()'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5
  assert a.mp == 5
  assert a.gallat is None
  assert a.output_map()['shrink'] == '5 5\n'
  assert repr(a) == 'Shrink(5)'
  assert repr(loads(dumps(a))) == repr(a)

  a.gallat = 10
  assert a.mp == 5
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '5 10\n'
  assert repr(a) == 'Shrink(5, 10)'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5, 
  assert all(array(a.mp) == [5])
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '0 10\n5 1 1\n'
  assert repr(a) == 'Shrink([5], 10)'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5, 6
  assert all(array(a.mp) == [5, 6])
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '0 10\n5 6 1\n'
  assert repr(a) == 'Shrink([5, 6], 10)'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5, 6, 7
  assert all(array(a.mp) == [5, 6, 7])
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '0 10\n5 6 7\n'
  assert repr(a) == 'Shrink([5, 6, 7], 10)'
  assert repr(loads(dumps(a))) == repr(a)


  a.raw = '5 5\n'
  assert a.mp == 5 and a.gallat is None
  a.raw = '5 10\n'
  assert a.mp == 5 and a.gallat == 10
  a.raw = '0 5\n5 1 1'
  assert all(array(a.mp) == [5]) and a.gallat is None
  a.raw = '0 10\n5 2 1'
  assert all(array(a.mp) == [5, 2]) and a.gallat == 10
  a.raw = '0 10\n5 2 10'
  assert all(array(a.mp) == [5, 2, 10]) and a.gallat == 10
  a.raw = '0 10\n1 1 1'
  assert all(array(a.mp) == 1) and a.gallat == 10
  assert repr(a) == 'Shrink(gallat=10)'

  try: a.gallat = 'a'
  except ValueError: pass
  else: raise Exception()
  try: a.mp = 'a'
  except ValueError: pass
  else: raise Exception()
  try: a.mp = range(4)
  except ValueError: pass
  else: raise Exception()

  assert a.output_map(structure=Molecule()) is None

def test_levshift():
  from pickle import loads, dumps
  from quantities import hartree, UnitQuantity, eV, kbar
  from lada.dftcrystal.electronic import LevShift, Electronic
  from lada.error import ValueError

  a = LevShift()
  assert a.output_map() is None
  assert a.units == UnitQuantity('decihartree', 0.1*hartree)
  assert a.lock is None and a.shift is None
  assert repr(a) == 'LevShift()'
  assert repr(loads(dumps(a))) == repr(a)

  a.shift = 0.005 * eV
  d = {'LevShift': LevShift}
  assert a.output_map() is None
  assert a.lock is None
  assert abs(a.shift - 0.005 * eV) < 1e-8
  assert abs(a.shift.magnitude - 0.005 * 1./0.1 * eV.rescale(hartree).magnitude) < 1e-8
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)
  
  a.shift = 2
  assert a.output_map() is None
  assert a.lock is None
  assert abs(a.shift - 0.2 * hartree) < 1e-8
  assert abs(a.shift.magnitude - 2) < 1e-8
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)

  a.lock = True
  assert a.raw == str(2) + ' ' + str(1)
  assert a.lock == True
  assert 'levshift' in a.output_map()
  assert a.output_map()['levshift'] == a.raw
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)

  a.lock = False
  assert a.raw == str(2) + ' ' + str(0)
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)

  a = Electronic()
  assert a.levshift.raw == ''
  a.levshift = 2, False
  assert a.levshift.raw == str(2) + ' ' + str(0)
  assert abs(a.levshift[0] - 0.2 * hartree) < 1e-8
  assert a.levshift[1] is False

  try: a.levshift = 1
  except ValueError: pass
  else: raise Exception()
  try: a.levshift = 1*kbar, False
  except: pass
  else: raise Exception()

if __name__ == '__main__': 
  test_shrink()
  test_levshift()
