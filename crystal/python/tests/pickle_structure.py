from pickle import loads, dumps
from numpy import all, abs, array
from lada.crystal import Structure
from lada.crystal._crystal import StructureStr, StructureVec, StructureSet

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell, weight=2, freeze=2, energy=-1e0, scale=2e0, name="2")
a.add_atom(0,0,0, "Si", site=0, freeze=2)\
          (0.25,0.25,0.25, "Ge", site=1, freeze=0)
a.magnetic = False
string = dumps(a)
b = loads(string)
assert type(b) == StructureStr
assert all(abs(a.cell-b.cell) < 1e-12)
assert abs(a.energy - b.energy) < 1e-12
assert abs(a.weight - b.weight) < 1e-12
assert abs(a.scale - b.scale) < 1e-12
assert a.freeze == b.freeze
assert a.name == b.name
assert hasattr(b, "magnetic")
assert a.magnetic == b.magnetic
assert all(abs(a[0].pos - b[0].pos) < 1e-12)
assert a[0].type == b[0].type
assert a[0].site == b[0].site
assert a[0].freeze == b[0].freeze
assert all(abs(a[1].pos - b[1].pos) < 1e-12)
assert a[1].type == b[1].type
assert a[1].site == b[1].site
assert a[1].freeze == b[1].freeze

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure( cell, weight=2, freeze=2, energy=-1e0,
               scale=2e0, name="2", kind="list" )
a.add_atom(0,0,0, "Si", "Ge", site=0, freeze=2)\
          (0.25,0.25,0.25, "Ge", site=1, freeze=0)
a.magnetic = False
string = dumps(a)
b = loads(string)
assert type(b) == StructureVec
assert all(abs(a.cell-b.cell) < 1e-12)
assert abs(a.energy - b.energy) < 1e-12
assert abs(a.weight - b.weight) < 1e-12
assert abs(a.scale - b.scale) < 1e-12
assert a.freeze == b.freeze
assert a.name == b.name
assert hasattr(b, "magnetic")
assert a.magnetic == b.magnetic
assert all(abs(a[0].pos - b[0].pos) < 1e-12)
assert a[0].type == b[0].type
assert a[0].site == b[0].site
assert a[0].freeze == b[0].freeze
assert all(abs(a[1].pos - b[1].pos) < 1e-12)
assert a[1].type == b[1].type
assert a[1].site == b[1].site
assert a[1].freeze == b[1].freeze

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure( cell, weight=2, freeze=2, energy=-1e0,
               scale=2e0, name="2", kind="set" )
a.add_atom(0,0,0, "Si", "Ge", site=0, freeze=2)\
          (0.25,0.25,0.25, "Ge", site=1, freeze=0)
a.magnetic = False
string = dumps(a)
b = loads(string)
assert type(b) == StructureSet
assert all(abs(a.cell-b.cell) < 1e-12)
assert abs(a.energy - b.energy) < 1e-12
assert abs(a.weight - b.weight) < 1e-12
assert abs(a.scale - b.scale) < 1e-12
assert a.freeze == b.freeze
assert a.name == b.name
assert hasattr(b, "magnetic")
assert a.magnetic == b.magnetic
assert all(abs(a[0].pos - b[0].pos) < 1e-12)
assert a[0].type == b[0].type
assert a[0].site == b[0].site
assert a[0].freeze == b[0].freeze
assert all(abs(a[1].pos - b[1].pos) < 1e-12)
assert a[1].type == b[1].type
assert a[1].site == b[1].site
assert a[1].freeze == b[1].freeze
