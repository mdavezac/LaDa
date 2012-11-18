""" Checks c++ to python and back.

    A cpp-created atom can be passed to python.
    A python-create atom can be passed to cpp and back.
    A python-subclassed atom can be passed to a cpp object and back. 
"""
def test(Class):
  import gc
  from numpy import all, abs
  from lada.crystal.cppwrappers import Atom
  from _atom_self import get_static_object, set_static_object, _atom
  
  # checks the static object is an Atom at start.
  a = get_static_object();
  assert a.__class__ is Atom
  assert a is _atom
  # checks the static object is always itself.
  b = get_static_object();
  assert a is b
  assert b is _atom
  # checks it can be changed and subclassed.
  a = Class(0.4,0.1,0.2, 'Au')
  set_static_object(a)
  c = get_static_object();
  assert c is not b
  assert c is a
  assert c.__class__ is Class
  assert all(abs(a.pos - [0.4, 0.1, 0.2]) < 1e-8) and a.type == 'Au'
  # checks same as above but with deletion.
  a = Class(0.4,0.1,0.2, 'Au', 'Pd', m=5)
  set_static_object(a)
  del a; del b; del c; gc.collect()
  c = get_static_object();
  assert c.__class__ is Class
  assert all(abs(c.pos - [0.4, 0.1, 0.2]) < 1e-8) and c.type == ['Au', 'Pd'] \
         and len(c.__dict__) == 1 and getattr(c, 'm', 0) == 5

if __name__ == "__main__": 
  from lada.crystal.cppwrappers import Atom
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  # tries to run test with normal class.
  test(Atom) 
  
  # tries to run test with other class. 
  # subclases an exception and makes sure it is called.
  # also checl passage through init.
  check_passage = False
  class MineErr(Exception): pass
  class Subclass(Atom):
    def __init__(self, *args, **kwargs):
      global check_passage
      check_passage = True
      super(Subclass, self).__init__(*args, **kwargs)
    def __call__(self): raise MineErr()

  test(Subclass)
  assert check_passage
 
  # checks that exception is called in __call__
  import gc
  from _atom_self import get_static_object, set_static_object
  gc.collect()
  try: get_static_object()()
  except MineErr: pass
  else: raise Exception()

  # tries to pass the wrong type object to the cpp wrapper.
  from lada.error import TypeError as LaDaTypeError
  class B(object): pass
  b = B()
  try: set_static_object(b)
  except LaDaTypeError: pass
  else: raise Exception()
