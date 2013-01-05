def test():
  from pylada.tools.input import Tree
  from pylada.error import ValueError

  a = Tree()
  assert len(a) == 0
  a = Tree(key=5)
  assert len(a) == 1
  assert a[0] == ('key', 5)
  assert a['key'] == 5
  for key in a.iterkeys(): assert key == 'key'
  for value in a.itervalues(): assert value == 5
  for key, value in a.iteritems(): assert key == 'key' and value == 5

  a = Tree(('key', 5), ('other', 'a'))
  assert len(a) == 2
  assert a[0] == ('key', 5)
  assert a['key'] == 5
  assert a[1] == ('other', 'a')
  assert a['other'] == 'a'
  iterator = a.iterkeys()
  assert iterator.next() == 'key'
  assert iterator.next() == 'other'
  try: iterator.next()
  except StopIteration: pass
  else: raise Exception()

  v = a.descend('branch', 'leaf')
  assert isinstance(v, Tree)
  assert isinstance(a['branch'], Tree)
  assert isinstance(a['branch'], Tree)
  assert isinstance(a['branch']['leaf'], Tree)
  assert a['branch']['leaf'] is v
  assert a[2][0] == 'branch'
  assert a[2][1] is a['branch']
  a['key'] = 6
  assert a['key'] == 6

  try: a[0] = 5
  except ValueError: pass
  else: raise Exception()
  try: a[0] = 5, 6
  except ValueError: pass
  else: raise Exception()


if __name__ == '__main__':
  test()
