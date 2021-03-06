""" Check memory deallocation.

    Creates cyclic references through an atom's dictionary.
    A number of instances are created and gets difference in memory hogged by python.
    These instances are then deleted and gc.collect is called.
    We then go through a loop where the instantiation/creation steps are
    perfomed again, with gc.collect to make sure garbage is deleted. 
    Finally, we check that memory is 10 times less than after first creation above.

    Note that python does not release memory it allocates. Hence, deletion and
    gc.collect must be inside the loop, or test could give false negative.

    The same is checked for a subclass of atom.
""" 
def test(Class): 
  import gc
  from os import getpid
  gc.set_debug(gc.DEBUG_OBJECTS | gc.DEBUG_UNCOLLECTABLE)

  id = getpid()
  def get_mem():
    from subprocess import Popen, PIPE
    output = Popen(["ps","-p", "{0}".format(id), '-o', 'rss'], stdout=PIPE).communicate()[0].splitlines()
    output = output[-1]
    return int(output)

  def mklist():
    result = [Class(0.1, 0.2, 0.5) for u in range(1000)]
    for b in result: b.m, b.type = ['Au', 'Pd']
    b = [(u.pos, u.type) for u in result]
    return result, b

  n = 10
  a = []
  startmem = get_mem()
  for i in range(n):
    a.append(mklist())
  mem = float(get_mem() - startmem) / float(n)
  assert mem > 0 # otherwise, test would be invalid.
  del a
  gc.collect()
   
  startmem = get_mem()
  for i in range(n*5): 
    a, b = mklist()
    # do deletion here, otherwise python might allocate extra memory to store our
    # objects, and the test would fail for reasons other than garbage collection.
    del a
    del b
    gc.collect()
  mem2 = float(get_mem() - startmem)
  assert mem2 < mem / 5.0
  assert len(gc.garbage) == 0

if __name__ == "__main__": 
  from pylada.crystal.cppwrappers import Atom
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  # tries to run test with normal class.
  test(Atom) 
 
  # tries to run test with other class. 
  # check passage through init.
  check_passage = False
  class Subclass(Atom):
    def __init__(self, *args, **kwargs):
      global check_passage
      check_passage = True
      super(Subclass, self).__init__(*args, **kwargs)

  test(Subclass)
  assert check_passage
