# Check memory is deallocated. 
from lada.crystal.cppwrappers import @PYTYPE@
import gc
from os import system, getpid
gc.set_debug(gc.DEBUG_OBJECTS | gc.DEBUG_UNCOLLECTABLE)

id = getpid()
def get_mem():
  from subprocess import Popen, PIPE
  output = Popen(["ps","--pid", "{0}".format(id), '-o', 'rss'], stdout=PIPE).communicate()[0].split('\n')[-2]
  return int(output)

def mklist():
  result = [@PYTYPE@() for u in range(1)]
  for b in result: b.m, b.type = [b], ['Au', 'Pd'] if "@PYTYPE@" == "AtomSequence" else 'Au'
  b = [(u.pos, u.type) for u in result]
  return result, b

n = 10
a = []
startmem = get_mem()
for i in range(n):
  a.append(mklist())
mem = float(get_mem() - startmem) / float(n)
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
assert mem2 < mem / 10.0
assert len(gc.garbage) == 0
