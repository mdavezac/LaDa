""" Checks atom methods and attributes. """
def test():
  """ Test atom initialization. """
  from _pyobject import represent, add_attribute, callme, equality

  assert represent("a") == repr("a");
  assert represent(1) == repr(1);

  class A(object):
    def __init__(self): 
      self.been_here = False
    def __call__(self): 
      self.been_here = True

  a = A()
  add_attribute(a, 'hello', 5)
  assert getattr(a, 'hello', 6) == 5
  assert not a.been_here
  callme(a)
  assert a.been_here
 
  assert equality(0, 0)
  assert not equality(1, 0)
  assert equality('1', '1')
  assert not equality('2', '1')

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  test()
