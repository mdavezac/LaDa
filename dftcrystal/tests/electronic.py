def test_shrink():
  from numpy import array, all
  from lada.dftcrystal.electronic import Shrink

  a = Shrink()
  assert a.mp == 1
  assert a.gallat is None
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'SHRINK'
  assert len(a.print_input().split('\n')[1].split()) == 2
  assert a.print_input().split('\n')[1].split()[0] == '1'
  assert a.print_input().split('\n')[1].split()[1] == '1'
  assert a.print_input().split('\n')[2] == ''
  assert repr(a) == 'Shrink()'

  a.mp = 5
  assert a.mp == 5
  assert a.gallat is None
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'SHRINK'
  assert len(a.print_input().split('\n')[1].split()) == 2
  assert a.print_input().split('\n')[1].split()[0] == '5'
  assert a.print_input().split('\n')[1].split()[1] == '5'
  assert a.print_input().split('\n')[2] == ''
  assert repr(a) == 'Shrink(5)'

  a.gallat = 10
  assert a.mp == 5
  assert a.gallat == 10
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'SHRINK'
  assert len(a.print_input().split('\n')[1].split()) == 2
  assert a.print_input().split('\n')[1].split()[0] == '5'
  assert a.print_input().split('\n')[1].split()[1] == '10'
  assert a.print_input().split('\n')[2] == ''
  assert repr(a) == 'Shrink(5, 10)'

  a.mp = 5, 
  assert all(array(a.mp) == [5])
  assert a.gallat == 10
  assert len(a.print_input().split('\n')) == 4
  assert a.print_input().split('\n')[0] == 'SHRINK'
  assert len(a.print_input().split('\n')[1].split()) == 2
  assert a.print_input().split('\n')[1].split()[0] == '0'
  assert a.print_input().split('\n')[1].split()[1] == '10'
  assert all(array(a.print_input().split('\n')[2].split(), dtype='int32') == [5, 1, 1])
  assert a.print_input().split('\n')[3] == ''
  assert repr(a) == 'Shrink([5], 10)'

  a.mp = 5, 6
  assert all(array(a.mp) == [5, 6])
  assert a.gallat == 10
  assert len(a.print_input().split('\n')) == 4
  assert a.print_input().split('\n')[0] == 'SHRINK'
  assert len(a.print_input().split('\n')[1].split()) == 2
  assert a.print_input().split('\n')[1].split()[0] == '0'
  assert a.print_input().split('\n')[1].split()[1] == '10'
  assert all(array(a.print_input().split('\n')[2].split(), dtype='int32') == [5, 6, 1])
  assert a.print_input().split('\n')[3] == ''
  assert repr(a) == 'Shrink([5, 6], 10)'

  a.mp = 5, 6, 7
  assert all(array(a.mp) == [5, 6, 7])
  assert a.gallat == 10
  assert len(a.print_input().split('\n')) == 4
  assert a.print_input().split('\n')[0] == 'SHRINK'
  assert len(a.print_input().split('\n')[1].split()) == 2
  assert a.print_input().split('\n')[1].split()[0] == '0'
  assert a.print_input().split('\n')[1].split()[1] == '10'
  assert all(array(a.print_input().split('\n')[2].split(), dtype='int32') == [5, 6, 7])
  assert a.print_input().split('\n')[3] == ''
  assert repr(a) == 'Shrink([5, 6, 7], 10)'


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

def test_levshift():
  from quantities import hartree, UnitQuantity, eV, kbar
  from lada.dftcrystal.electronic import LevShift, Electronic
  from lada.error import ValueError

  a = LevShift()
  assert a.print_input() is None
  assert a.units == UnitQuantity('decihartree', 0.1*hartree)
  assert a.lock is None and a.shift is None

  a.shift = 0.005 * eV
  assert a.print_input() is None
  assert a.lock is None
  assert abs(a.shift - 0.005 * eV) < 1e-8
  assert abs(a.shift.magnitude - 0.005 * 1./0.1 * eV.rescale(hartree).magnitude) < 1e-8
  
  a.shift = 2
  assert a.print_input() is None
  assert a.lock is None
  assert abs(a.shift - 0.2 * hartree) < 1e-8
  assert abs(a.shift.magnitude - 2) < 1e-8

  a.lock = True
  assert a.raw == str(float(2)) + ' ' + str(1)
  assert a.lock == True
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'LEVSHIFT'
  assert a.print_input().split('\n')[1] ==  a.raw 
  assert a.print_input().split('\n')[-1] == ''

  a.lock = False
  assert a.raw == str(float(2)) + ' ' + str(0)

  a = Electronic()
  assert a.levshift.raw == ''
  a.levshift = 2, False
  assert a.levshift.raw == str(float(2)) + ' ' + str(0)
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