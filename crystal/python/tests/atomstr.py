import gc
from numpy import all
from lada.crystal.cppwrappers.atom import AtomStr

# # Try correct initialization. Check for garbage collection.
# a = AtomStr()
# assert a.type == "" and a.site == -1 and a.freeze == 0 and all(a.pos == 0) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)
# else: print "Garbage collection disabled."

# a = AtomStr('Au')
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# a = AtomStr(0.1, 0.1, 0.1, 'Au')
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# a = AtomStr([0.1, 0.1, 0.1], 'Au')
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# a = AtomStr([0.1, 0.1, 0.1], type='Au')
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# a = AtomStr(type='Au', pos=[0.1, 0.1, 0.1])
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# a = AtomStr('Au', pos=[0.1, 0.1, 0.1])
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# a = AtomStr('Au', pos=[0.1, 0.1, 0.1], m=1.3)
# assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1)\
#        and len(a.__dict__) == 1 and abs(a.__dict__.get('m', 0) - 1.3) < 1e-12
# if gc.isenabled():
#   assert gc.is_tracked(a)
#   del a
#   gc.collect()
#   assert gc.get_count() == (0, 0, 0)

# # Try incorrect initialization. Check for garbage collection.
# try: a = AtomStr([0.1, 0.1])
# except TypeError: pass
# else: raise RuntimeError("Should have failed.")
# finally:
#   if gc.isenabled():
#     gc.collect()
#     assert gc.get_count() == (0, 0, 0)

# try: a = AtomStr([0.1, 0.1, 0.1, 0.1])
# except TypeError: pass
# else: raise RuntimeError("Should have failed.")
# finally:
#   if gc.isenabled():
#     gc.collect()
#     assert gc.get_count() == (0, 0, 0)

# try: a = AtomStr(0.1, 0.1, 0.1, 0.1)
# except TypeError: pass
# else: raise RuntimeError("Should have failed.")
# finally:
#   if gc.isenabled():
#     gc.collect()
#     assert gc.get_count() == (0, 0, 0)

# try: a = AtomStr(0.1, 'Au', 0.1, 0.1)
# except TypeError: pass
# else: raise RuntimeError("Should have failed.")
# finally:
#   if gc.isenabled():
#     gc.collect()
#     assert gc.get_count() == (0, 0, 0)

# try: a = AtomStr('Au', type='Au')
# except TypeError: pass
# else: raise RuntimeError("Should have failed.")
# finally:
#   if gc.isenabled():
#     gc.collect()
#     assert gc.get_count() == (0, 0, 0)

# try: a = AtomStr(0.1, 0.1, 0.1, pos=[0.1, 0.1, 0.1])
# except TypeError: pass
# else: raise RuntimeError("Should have failed.")
# finally:
#   if gc.isenabled():
#     gc.collect()
#     assert gc.get_count() == (0, 0, 0)

from time import sleep
for i in range(10):
  a = [AtomStr() for u in range(10)]
  del a
  sleep(2)

