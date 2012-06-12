
def test():
  from lada.functools.makeclass import makeclass
  import dummy

  base = dummy.Base(1, 'b')

  Functional = makeclass('Functional', dummy.Base, dummy.iter_func, call=dummy.func)
  functional = Functional(b=-5, other=False)
  assert getattr(functional, 'other', True) == False
  assert functional(True) == False
  assert functional.a == 2
  assert dummy.iterator == 4
  assert functional.b == -5

  functional = Functional(b=-5, other=False, copy=base)
  assert functional.a == 1
  assert getattr(functional, 'other', True) == False
  assert functional(True) == False
  assert functional.a == 1
  assert dummy.iterator == 4
  assert functional.b == -5

  Functional = makeclass('Functional', dummy.Base, dummy.iter_func)
  functional = Functional(b=-5, other=False)
  assert getattr(functional, 'other', True) == False
  assert functional(True) is None
  assert functional.a == 2
  assert dummy.iterator == 4
  assert functional.b == -5

  functional = Functional(b=-5, other=False, copy=base)
  assert functional.a == 1
  assert getattr(functional, 'other', True) == False
  assert functional(True) is None
  assert functional.a == 1
  assert dummy.iterator == 4
  assert functional.b == -5

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  test()
