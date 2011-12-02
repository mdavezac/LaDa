from lada.crystal.cppwrappers.atom import AtomSequence
import gc
from os import system, getpid

print "HERE"
a = AtomSequence()
print "THERE"
del a
gc.collect()
print "HERE"
 
id = getpid()
for i in range(100): 
  a = [AtomSequence() for u in range(100000)]
  for b in a: b.m = 0
  del a
  del b
  gc.collect()
  system("ps --pid {0} u".format(id))
# print h.heap()
