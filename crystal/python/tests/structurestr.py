import gc
from numpy import all, array
from lada.crystal.cppwrappers import StructureStr

# Try correct initialization.
a = StructureStr()
assert all(abs(a.cell - [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "" and abs(a.scale - 1e0) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
a.cell = [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "" and abs(a.scale - 1e0) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
for i in range(3):
  for j in range(3):
    a.cell[i, j] = 5
    assert abs(a.cell[i, j] - 5) < 1e-12
a.cell = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
assert all(abs(a.cell - [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "" and abs(a.scale - 1e0) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
a = StructureStr([[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]])
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "" and abs(a.scale - 1e0) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
a = StructureStr(cell=[[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]])
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "" and abs(a.scale - 1e0) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
a = StructureStr(0.5, cell=[[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]])
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "" and abs(a.scale - 0.5) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
a = StructureStr(0.5, "what", cell=[[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]])
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 0 and a.name == "what" and abs(a.scale - 0.5) < 1e-12   \
       and abs(a.energy - 0e0) < 1e-12 and abs(a.weight - 1e0) < 1e-12
a = StructureStr(0.5, "what", cell=[[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]], freeze=1, weight=0.5, energy=0.2)
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 1 and a.name == "what" and abs(a.scale - 0.5) < 1e-12   \
       and abs(a.energy - 0.2) < 1e-12 and abs(a.weight - 0.5) < 1e-12
a = StructureStr(0.5, "what", cell=[[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]], freeze=1, weight=0.5, energy=0.2, m=5)
assert all(abs(a.cell - [[1, 1, 0], [0, 1, 0.4], [0, 0.1, 1]]) < 1e-12)        \
       and a.freeze == 1 and a.name == "what" and abs(a.scale - 0.5) < 1e-12   \
       and abs(a.energy - 0.2) < 1e-12 and abs(a.weight - 0.5) < 1e-12         \
       and len(a.__dict__) == 1 and getattr(a, 'm', 0) == 5
assert repr(a) == "structure = StructureStr(0.5, 'what', weight=0.5, energy=0.2, freeze=1, m=5)\n"\
                  "structure.cell = [ [1, 0, 0],\\\n"\
                  "                   [1, 1, 0.1],\\\n"\
                  "                   [0, 0.4, 1] ]"
class B(StructureStr):
  def __new__(self, *args, **kwargs):
    print "THERE"
    return super(B, self).__new__(self, *args, **kwargs)
  def __init__(self): 
    print "HERE"
    super(B, self).__init__()
  def __getitem__(self, *args, **kwargs):
    print ">> ", args
    print kwargs
print B.__name__
print B.__bases__
a = B()
print a.__class__.__name__
print len(a)

