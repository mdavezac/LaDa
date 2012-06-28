def test_boolkeyword():
  from lada.dftcrystal.input import BoolKeyword

  class Dummy(BoolKeyword): 
    keyword = 'whatever'
    pass

  a = Dummy('whatever')
  a.value = None
  assert a.value is False
  assert a.print_input() is None
  a.value = 1
  assert a.value is True
  assert a.print_input() == 'WHATEVER\n'

def test_valuekeyword():
  from lada.dftcrystal.input import ValueKeyword

  a = ValueKeyword('whatever')
  assert a.keyword == 'whatever'
  assert a.value is None
  assert a.print_input() is None

  a.raw = '5'
  assert type(a.value) is int and a.value == 5 and a.raw == '5'
  a.raw = '5.0'
  assert type(a.value) is float and abs(a.value - 5.0) < 1e-8 and a.raw == str(5.0)
  a.raw = 'five'
  assert type(a.value) is str and a.value == 'five' and a.raw == 'FIVE'
  a.raw = '5 5.0 five'
  assert hasattr(a.value, '__iter__') and hasattr(a.value, '__len__')
  assert len(a.value) == 3
  assert type(a.value[0]) is int and a.value[0] == 5 
  assert type(a.value[1]) is float and abs(a.value[1] - 5.0) < 1e-8
  assert type(a.value[2]) is str and a.value[2] == 'five'
  assert type(a.raw) is str and len(a.raw.split('\n')) == 1 and len(a.raw.split()) == 3
  assert a.raw.split()[0] == '5' and a.raw.split()[1] == str(5.0) and a.raw.split()[2] == 'FIVE' 

  string = a.print_input()
  assert len(string.split('\n')) == 3
  line0, line1, line2 = string.split('\n')
  assert line2 == ''
  assert line0 == 'WHATEVER'
  line1 = line1.split()
  assert len(line1) == 3
  assert line1[0] == '5' and line1[1] == str(5.0) and line1[2] == 'FIVE' 
  
  a.value = None
  assert a.value is None
  assert a.print_input() is None

def test_typedkeyword():
  from lada.dftcrystal.input import TypedKeyword

  a = TypedKeyword('whatever', int)
  assert a.keyword == 'whatever'
  assert a.value is None
  assert a.print_input() is None
  a.raw = '5'
  assert type(a.value) is int and a.value == 5 and a.raw == '5'
  a.value = '5'
  assert type(a.value) is int and a.value == 5 and a.raw == '5'
  a.value = None
  assert a.value is None
  assert a.print_input() is None
  a.value = 5.1
  assert type(a.value) is int and a.value == 5 and a.raw == '5'
  assert len(a.print_input().split('\n')) == 3 
  assert a.print_input().split('\n')[0] == 'WHATEVER' 
  assert a.print_input().split('\n')[1] == '5'
  assert a.print_input().split('\n')[2] == ''
  try: a.value = 'five'
  except: pass
  else: raise Exception()
  
  a = TypedKeyword('whatever', float)
  assert a.keyword == 'whatever'
  assert a.value is None
  assert a.print_input() is None
  a.raw = '5.1'
  assert type(a.value) is float and abs(a.value - 5.1) < 1e-8 and a.raw == str(5.1)
  a.value = '5.1'
  assert type(a.value) is float and abs(a.value - 5.1) < 1e-8 and a.raw == str(5.1)
  a.value = None
  assert a.value is None
  assert a.print_input() is None
  a.value = 5.1
  assert type(a.value) is float and abs(a.value - 5.1) < 1e-8 and a.raw == str(5.1)
  assert len(a.print_input().split('\n')) == 3 
  assert a.print_input().split('\n')[0] == 'WHATEVER' 
  assert a.print_input().split('\n')[1] == str(5.1)
  try: a.value = 'five'
  except: pass
  else: raise Exception()
  
  a = TypedKeyword('whatever', str)
  assert a.keyword == 'whatever'
  assert a.value is None
  assert a.print_input() is None
  a.raw = '5.1'
  assert type(a.value) is str and a.value == '5.1' and a.raw == '5.1'
  a.value = '5.1'
  assert type(a.value) is str and a.value == '5.1' and a.raw == '5.1'
  a.value = None
  assert a.value is None
  assert a.print_input() is None
  a.value = 5.1
  assert type(a.value) is str and a.value == str(5.1) and a.raw == str(5.1)
  assert len(a.print_input().split('\n')) == 3 
  assert a.print_input().split('\n')[0] == 'WHATEVER' 
  assert a.print_input().split('\n')[1] == str(5.1)
  a.value = 'five'
  assert type(a.value) is str and a.value == 'five' and a.raw == 'FIVE'

  a = TypedKeyword('whatever', [int])
  a.raw = '5 5'
  assert hasattr(a.value, '__iter__') 
  assert len(a.value) == 2
  assert all(type(v) is int for v in a.value)
  assert all(v == 5 for v in a.value)
  a.value = 5, 5.3, 5.2
  assert hasattr(a.value, '__iter__') 
  assert len(a.value) == 3
  assert all(type(v) is int for v in a.value)
  assert all(v == 5 for v in a.value)
  assert all(v == '5' for v in a.raw.split())
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'WHATEVER' 
  assert all(v == '5' for v in a.print_input().split('\n')[1].split())
  assert a.print_input().split('\n')[2] == '' 
  a.raw = ''
  assert a.value is None
  a.value = []
  assert a.value is None
  try: a.raw = '0 FIVE'
  except: pass
  else: raise Exception()
  
  a = TypedKeyword('whatever', [str, int])
  a.raw = 'five 5'
  assert hasattr(a.value, '__iter__') 
  assert len(a.value) == 2
  assert type(a.value[0]) is str and a.value[0] == 'five'
  assert type(a.value[1]) is int and a.value[1] == 5
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'WHATEVER' 
  assert a.print_input().split('\n')[1].split()[0] == 'FIVE' 
  assert a.print_input().split('\n')[1].split()[1] == '5' 
  assert a.print_input().split('\n')[2] == '' 
  try: a.value = 'five', 6, 7
  except: pass
  else: raise Exception()
  try: a.value = 'five 5'
  except: pass
  else: raise Exception()
  try: a.raw = '0 five'
  except: pass
  else: raise Exception()

def test_choicekeyword():
  from lada.dftcrystal.input import Choice

  a = Choice(['five', 5], keyword='whatever')
  assert a.value is None
  assert a.print_input() is None

  a.value = 5
  assert a.value == 5
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'WHATEVER'
  assert a.print_input().split('\n')[1] == '5'
  assert a.print_input().split('\n')[2] == ''

  a.value = 'five'
  assert a.value == 'five'
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'WHATEVER'
  assert a.print_input().split('\n')[1] == 'FIVE'
  assert a.print_input().split('\n')[2] == ''
  a.value = 'FIVE'
  assert a.value == 'five'
  assert len(a.print_input().split('\n')) == 3
  assert a.print_input().split('\n')[0] == 'WHATEVER'
  assert a.print_input().split('\n')[1] == 'FIVE'
  assert a.print_input().split('\n')[2] == ''

  try: a.value = 6
  except: pass
  else: raise Exception()
  try: a.value = 'six'
  except: pass
  else: raise Exception()
 
if __name__ == "__main__":
  test_choicekeyword()
  test_valuekeyword()
  test_typedkeyword()
  test_boolkeyword()

