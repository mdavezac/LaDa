from lada.crystal._crystal import AtomStr, AtomVec, AtomSet, FreezeAtom, SpecieSet, VectorStr
from lada import error
from numpy import all, abs

a = AtomSet()
assert all(abs(a.pos) < 1e-12)
assert len(a.type) == 0
assert type(a.type) == SpecieSet
assert a.freeze == 0
assert a.site == -1
a = AtomSet(0, 1, 2, "Au", "Pd", site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == ["Au", "Pd"]
assert type(a.type) == SpecieSet
assert a.freeze == FreezeAtom.x
assert a.site == 0
a = AtomSet(0, 1, 2, SpecieSet(["Au", "Pd"]), site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == ["Au", "Pd"]
assert type(a.type) == SpecieSet
assert a.freeze == FreezeAtom.x
assert a.site == 0
a = AtomSet(position=[0, 1, 2], type=["Au", "Pd"], site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == ["Au", "Pd"]
assert type(a.type) == SpecieSet
assert a.freeze == FreezeAtom.x
assert a.site == 0
a.pos[0] = -1.5
assert abs(a.pos[0] + 1.5) < 1e-12 
a.pos[:] = -2.2
assert all(abs(a.pos + 2.2)) < 1e-12 

a = AtomVec()
assert all(abs(a.pos) < 1e-12)
assert len(a.type) == 0
assert type(a.type) == VectorStr
assert a.freeze == 0
assert a.site == -1
a = AtomVec(0, 1, 2, "Au", "Pd", site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == ["Au", "Pd"]
assert type(a.type) == VectorStr
assert a.freeze == FreezeAtom.x
assert a.site == 0
a = AtomVec(0, 1, 2, ["Au", "Pd"], site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == ["Au", "Pd"]
assert type(a.type) == VectorStr
assert a.freeze == FreezeAtom.x
assert a.site == 0
a = AtomVec(position=[0, 1, 2], type=["Au", "Pd"], site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == ["Au", "Pd"]
assert type(a.type) == VectorStr
assert a.freeze == FreezeAtom.x
assert a.site == 0
a.pos[0] = -1.5
assert abs(a.pos[0] + 1.5) < 1e-12 
a.pos[:] = -2.2
assert all(abs(a.pos + 2.2)) < 1e-12 

a = AtomStr()
assert all(abs(a.pos) < 1e-12)
assert len(a.type) == 0
assert type(a.type) == str
assert a.freeze == 0
assert a.site == -1
a = AtomStr(0, 1, 2, "Au", site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == "Au"
assert type(a.type) == str
assert a.freeze == FreezeAtom.x
assert a.site == 0
a = AtomStr(position=[0, 1, 2], type="Au", site=0, freeze=1)
assert abs(a.pos[0]) < 1e-12
assert abs(a.pos[1] - 1) < 1e-12
assert abs(a.pos[2] - 2) < 1e-12
assert a.type == "Au"
assert type(a.type) == str
assert a.freeze == FreezeAtom.x
assert a.site == 0
a.pos[0] = -1.5
assert abs(a.pos[0] + 1.5) < 1e-12 
a.pos[:] = -2.2
assert all(abs(a.pos + 2.2)) < 1e-12 
try: a.type = ["Au", "Pd"]
except error.ValueError: pass
try: a.type = ["Au", "Pd"]
except ValueError: pass
