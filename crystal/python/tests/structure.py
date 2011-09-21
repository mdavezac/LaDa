from numpy import all, abs, array, any
from lada.crystal import Structure
from lada.crystal._crystal import StructureStr, StructureVec, StructureSet

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell)
assert type(a) == StructureStr
assert a.kind == "scalar"
assert all(abs(cell - a.cell) < 1e-12)
cell[:] = 0e0
assert any(abs(cell - a.cell) > 1e-12)
assert abs(a.scale - 1e0)  < 1e-12
assert abs(a.weight - 1e0) < 1e-12
assert abs(a.energy - 0e0) < 1e-12
assert len(a) == 0
a.add_atom(0,0,0, "Si")\
          (0.25,0.25,0.25, "Ge")
assert len(a) == 2
assert all(abs(a[0].pos) < 1e-12)
assert a[0].site == -1
assert a[0].freeze == 0
assert a[0].type == "Si"
assert all(abs(a[1].pos - 0.25) < 1e-12)
assert a[1].site == -1
assert a[1].freeze == 0
assert a[1].type == "Ge"

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell, kind="list", weight=0e0)
assert type(a) == StructureVec
assert a.kind == "list"
assert all(abs(cell - a.cell) < 1e-12)
cell[:] = 0e0
assert any(abs(cell - a.cell) > 1e-12)
assert abs(a.scale - 1e0)  < 1e-12
assert abs(a.weight - 0e0) < 1e-12
assert abs(a.energy - 0e0) < 1e-12
assert len(a) == 0
a.add_atom(0,0,0, "Si")\
          (0.25,0.25,0.25, "Si", "Ge")
assert len(a) == 2
assert all(abs(a[0].pos) < 1e-12)
assert a[0].site == -1
assert a[0].freeze == 0
assert a[0].type == ["Si"]
assert all(abs(a[1].pos - 0.25) < 1e-12)
assert a[1].site == -1
assert a[1].freeze == 0
assert a[1].type == ["Si", "Ge"]

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell, kind="set", weight=0e0)
assert type(a) == StructureSet
assert a.kind == "set"
assert all(abs(cell - a.cell) < 1e-12)
cell[:] = 0e0
assert any(abs(cell - a.cell) > 1e-12)
assert abs(a.scale - 1e0)  < 1e-12
assert abs(a.weight - 0e0) < 1e-12
assert abs(a.energy - 0e0) < 1e-12
assert len(a) == 0
a.add_atom(0,0,0, "Si")\
          (0.25,0.25,0.25, "Si", "Ge", "Si")
assert len(a) == 2
assert all(abs(a[0].pos) < 1e-12)
assert a[0].site == -1
assert a[0].freeze == 0
assert a[0].type == ["Si"]
assert all(abs(a[1].pos - 0.25) < 1e-12)
assert a[1].site == -1
assert a[1].freeze == 0
assert a[1].type == ["Si", "Ge"]
