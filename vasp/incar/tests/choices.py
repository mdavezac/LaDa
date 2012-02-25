def test():
  from pickle import loads, dumps
  from lada.vasp.incar._params import Choices

  a = Choices('algo', {'A': ['aa', 0], 'B': ['bb', 1]})
  a.value = 'a'
  assert a.value == 'A'
  assert loads(dumps(a)).incar_string() == 'ALGO = A'
  a.value = 'aa'
  assert a.value == 'A'
  assert loads(dumps(a)).incar_string() == 'ALGO = A'
  a.value = 0
  assert a.value == 'A'
  assert loads(dumps(a)).incar_string() == 'ALGO = A'
  a.value = 'b'
  assert a.value == 'B'
  assert loads(dumps(a)).incar_string() == 'ALGO = B'
  a.value = 'bb'
  assert a.value == 'B'
  assert loads(dumps(a)).incar_string() == 'ALGO = B'
  a.value = 1
  assert a.value == 'B'
  assert loads(dumps(a)).incar_string() == 'ALGO = B'
  try: a.value = 2
  except: pass
  else: raise RuntimeError()
  a.value = None
  assert a.incar_string() is None
  assert repr(a) == "Choices('algo', {'A': ['aa', 0, 'a'], 'B': ['bb', 1, 'b']}, None)"

  a = Choices('algo', {'A': ['aa', 0], 'B': ['bb', 1]}, 'b')
  assert a.value == 'B'
  a = Choices('algo', {'A': ['aa', 0], 'B': ['bb', 1]}, 'bb')
  assert a.value == 'B'
  a = Choices('algo', {'A': ['aa', 0], 'B': ['bb', 1]}, 1)
  assert a.value == 'B'
  try: a = Choices('algo', {'A': ['aa', 0], 'B': ['bb', 1]}, 2)
  except: pass
  else: raise RuntimeError()

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

