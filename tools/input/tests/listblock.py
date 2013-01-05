def test():
  from pickle import loads, dumps
  from pylada.tools.input import ListBlock, BaseKeyword

  a = ListBlock()
  d = {'ListBlock': a.__class__}
  assert len(a) == 0
  assert isinstance(eval(repr(a), d), ListBlock)
  assert len(eval(repr(a), d)) == 0
  assert isinstance(loads(dumps(a)), ListBlock)
  assert len(loads(dumps(a))) == 0
  a.append('hello', True)
  assert len(a) == 1
  assert a[0].__class__ is BaseKeyword
  assert a[0].keyword == 'hello'
  assert a[0].raw == True
  assert a.output_map() == [('hello', True)]
  assert isinstance(eval(repr(a), d), ListBlock)
  assert len(eval(repr(a), d)) == 1
  assert eval(repr(a), d)[0].__class__ is BaseKeyword
  assert eval(repr(a), d)[0].keyword == 'hello'
  assert eval(repr(a), d)[0].raw == True
  assert len(loads(dumps(a))) == 1
  assert loads(dumps(a))[0].__class__ is BaseKeyword
  assert loads(dumps(a))[0].keyword == 'hello'
  assert loads(dumps(a))[0].raw == True
  b = ListBlock()
  b.read_input(a.output_map())
  assert repr(b) == repr(a)
  

if __name__ == '__main__':
  test()
