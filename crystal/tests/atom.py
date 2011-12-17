""" Checks atom methods and attributes. """
def test_init(Class):
  """ Test atom initialization. """
  from numpy import all, abs, array

  # Try correct initialization. Check for garbage collection.
  a = Class()
  assert a.type is None and all(abs(a.pos) < 1e-8) and len(a.__dict__) == 0
  assert a.type is None and all(abs(a.pos) < 1e-8) and len(a.__dict__) == 0
  a = Class(0.1, 0.1, 0.1, 'Au')
  assert a.type == "Au" and all(abs(a.pos-0.1) < 1e-8) and len(a.__dict__) == 0
  a = Class(0.1, 0.1, 0.1, type=['Au', 'Pd'])
  assert a.type == ["Au", 'Pd'] and all(abs(a.pos-0.1) < 1e-8) and len(a.__dict__) == 0
  a = Class(type='Au', pos=[0.1, 0.1, 0.1])
  assert a.type == "Au" and all(abs(a.pos-0.1) < 1e-8) and len(a.__dict__) == 0
  a = Class(type='Au', pos=[0.1, 0.1, 0.1], m=5)
  assert a.type == "Au" and all(abs(a.pos-0.1) < 1e-8) and len(a.__dict__) == 1 and getattr(a, 'm', 3) == 5
  a = Class(0.1, 0.1, 0.1, 0.1, 0.1)
  assert all(abs(array(a.type) - 0.1) < 1e-8) and all(abs(a.pos-0.1) < 1e-8) and len(a.__dict__) == 0
  l = [None]
  a = Class(0.1, 0.1, 0.1, l)
  assert a.type is l
  assert all(abs(a.pos-0.1) < 1e-8) and len(a.__dict__) == 0
  a.pos[0] = 0.2
  a.pos[1] = 0.3
  a.pos[2] = 0.4
  assert all(abs(a.pos-[0.2, 0.3, 0.4]) < 1e-8)
  
def test_fail_init(Class):
  """ Test failures during initialization. """
  # Try incorrect initialization
  try: a = Class(0.1, 0.1)
  except TypeError: pass
  else: raise RuntimeError("Should have failed.")
  try: a = Class(0.1, 'Au', 0.1, 0.1)
  except TypeError: pass
  else: raise RuntimeError("Should have failed.")
  try: a = Class(0.1, 0.1, 0.1, 'Au', type='Au')
  except TypeError: pass
  else: raise RuntimeError("Should have failed.")
  try: a = Class(0.1, 0.1, 0.1, pos=[0.1, 0.1, 0.1])
  except TypeError: pass
  else: raise RuntimeError("Should have failed.")

def test_repr(Class):
  """ Test representability. """
  assert repr(Class(type='Au', pos=[1, 1, 1], m=1)) == "{0.__name__}(1, 1, 1, 'Au', m=1)".format(Class)
  assert str(Class(type='Au', pos=[1, 1, 1], site=1)) == "{0.__name__}(1, 1, 1, 'Au', site=1)".format(Class)

def test_copy(Class):
  """ Test copy and deep copy. """
  from copy import copy, deepcopy
  b = Class(0,0,0, 'Au', m=0)
  a = copy(b)
  b.type = 'Pd'
  assert a is b
  a = deepcopy(b)
  b.type = 'Pd'
  b.pos += 1
  del b.m 
  assert a.type == "Pd" and all(abs(a.pos-0.0) < 1e-8) and len(a.__dict__) == 1 and getattr(a, 'm', 1) == 0
  assert a.__class__ is Class

def test_todict(Class):
  """ Test to_dict member. """
  a = Class(0,0,0, 'Au', m=0)
  a = a.to_dict()
  assert a['type'] == "Au" and all(a['pos']== 0) and a['m'] == 0

def test_pickle(Class):
  """ Test pickling. """
  from pickle import loads, dumps
  a = Class(pos=[0, 1, 2], type="Au", m=6)
  b = loads(dumps(a))
  assert b.__class__ is a.__class__
  assert all(abs(b.pos - [0, 1, 2]) < 1e-12) and b.type == 'Au' \
         and len(b.__dict__) == 1 and b.__dict__.get('m', 0) == 6
  b = loads(dumps((a, a)))
  assert b[0] is b[1]

if __name__ == "__main__":
  from lada.crystal.cppwrappers import Atom
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  # tries to run test with normal class.
  test_init(Atom) 
  test_fail_init(Atom) 
  test_repr(Atom) 
  test_copy(Atom) 
  test_todict(Atom) 
  test_pickle(Atom) 
 
  # tries to run test with other class. 
  # check passage through init.
  check_passage = False
  class Subclass(Atom):
    def __init__(self, *args, **kwargs):
      global check_passage
      check_passage = True
      super(Subclass, self).__init__(*args, **kwargs)

  test_init(Subclass) 
  test_fail_init(Subclass) 
  test_repr(Subclass) 
  test_copy(Subclass) 
  test_todict(Subclass) 
  test_pickle(Subclass) 
  assert check_passage

