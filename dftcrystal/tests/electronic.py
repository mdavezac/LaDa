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

if __name__ == '__main__': 
  test_shrink()
