def test():
  from lada.jobfolder.forwarding_dict import ForwardingDict
  from collections import namedtuple

  class A(object):
    def __init__(self, this, that):
      self.this = this
      self.that = that
    def __eq__(self, b):
      return b.__class__ is self.__class__ and b.this == self.this and b.that == self.that
    def __neq__(self, b): return not self.__eq__(b)
    def __repr__(self):  return "A({0.this}, {0.that})".format(self)
  first  = A(0, A(5, A('a', 'b')))
  second = A(1, A(6, A('c', 'd')))
  third  = A(0, A(5, A('a', 'b')))
  assert first == third
  third.that.that.that = 'd'
  assert first != third

  d = ForwardingDict(ordered=True, readonly=True)
  d['first'] = first
  d['second'] = second
  assert     d['first'].this == 0 and d['first'].that.this == 5\
         and d['first'].that.that.this == 'a' and d['first'].that.that.that == 'b'
  assert     d['second'].this == 1 and d['second'].that.this == 6\
         and d['second'].that.that.this == 'c' and d['second'].that.that.that == 'd'
  assert repr(d) == "{\n  'second': A(1, A(6, A(c, d))),\n  'first':  A(0, A(5, A(a, b))),\n}"
  assert repr(d.this) == "{\n  'second': 1,\n  'first':  0,\n}"
  assert repr(d.that) == "{\n  'second': A(6, A(c, d)),\n  'first':  A(5, A(a, b)),\n}"
  assert repr(d.that.this) == "{\n  'second': 6,\n  'first':  5,\n}"
  assert repr(d.that.that) == "{\n  'second': A(c, d),\n  'first':  A(a, b),\n}"
  assert repr(d.that.that.this) == "{\n  'second': 'c',\n  'first':  'a',\n}"
  assert repr(d.that.that.that) == "{\n  'second': 'd',\n  'first':  'b',\n}"

  for key, value in d.iteritems():
    assert {'first': first, 'second': second }[key] == value
  for key, value in d.this.iteritems():
    assert {'first': first.this, 'second': second.this }[key] == value
  for key, value in d.that.iteritems():
    assert {'first': first.that, 'second': second.that }[key] == value
  for key, value in d.that.this.iteritems():
    assert {'first': first.that.this, 'second': second.that.this }[key] == value
  for key, value in d.that.that.iteritems():
    assert {'first': first.that.that, 'second': second.that.that }[key] == value
  for key, value in d.that.that.this.iteritems():
    assert {'first': first.that.that.this, 'second': second.that.that.this }[key] == value
  for key, value in d.that.that.that.iteritems():
    assert {'first': first.that.that.that, 'second': second.that.that.that }[key] == value

  try: d.this.that
  except: pass
  else: raise RuntimeError()
  try: d.this = 8
  except: pass
  else: raise RuntimeError()
  try: del d.this
  except: pass
  else: raise RuntimeError()

  d.readonly = False
  assert all([u != 8 for u in d.this.itervalues()])
  d.this = 8
  assert all([u == 8 for u in d.this.itervalues()])
  assert all([u != 8 for u in d.that.this.itervalues()])
  d.that.this = 8
  assert all([u == 8 for u in d.that.this.itervalues()])
  del d.that.this
  assert not hasattr(first.that, 'this')
  assert not hasattr(second.that, 'this')
  assert hasattr(first.that.that, 'this')  and hasattr(first, 'this') 
  assert hasattr(second.that.that, 'this') and hasattr(second, 'this')
  try: print d.that.this
  except: pass
  else: raise RuntimeError()

  d.naked_end = False
  first.that.this = 5
  assert d.that.this.values()[0] == first.that.this
  d.naked_end = True
  assert d.that.this == first.that.this
  
  d.only_existing = True
  try: d.that.other = True
  except: pass
  else: raise RuntimeError()
  d.only_existing = False
  d.that.other = True
  assert getattr(first.that, 'other', False) == True
  assert getattr(second.that, 'other', False) == True

  del first.that.that
  d.that.that.other = True
  assert getattr(second.that.that, 'other', False) == True


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
